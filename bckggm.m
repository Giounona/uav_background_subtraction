function [model,idx,res,Mymodel] = bckggm(im,model,idx,cfg,rot_ang,dx,dy,Mymodel,r_h, mask)
% BCKGGM Adaptive background modeling by using a mixture of Gaussians
% CMP Vision Algorithms http://visionbook.felk.cvut.cz
% Tomas Svoboda, 2007
% 
%

%Modified by Giounona Tzanidou for background modelling for aerial video
%surveillance
%date: 3-07-2015



% This function computes a one-frame update
% of the background model and performs motion detection; it is an
% implementation of . Each pixel is modeled by
% a weighted mixture of Gaussians whose parameters and weights are
% continuously updated throughout the image sequence.
onemat=ones(cfg.K,1);
%Newmodel=model;%%*******************************************************
NewMymodel=Mymodel;
im = double(im);

% If there is no model, initialize it. See model_init.
if isempty(Mymodel)
  [model,idx,Mymodel] = model_init2( im, cfg );
end

Him=size(im,1);
Wim=size(im,2);
halfH=(Him+1)/2;
halfW=(Wim+1)/2;


% Memory allocation for the results
res.segm = zeros( [size(im,1),size(im,2)], 'uint8' );
res.bckg = zeros( size(im) );
layers = size( im, 3 );
%siz = size( model(1,1,:,:) );
% Inspect each pixel and update its background model
for i = 1:size(im,1)
  for j = 1:size(im,2)
% Prepare the data.
% x is a vector containing layer values, for an RGB image thus x=[r,g,b].
% m
% is a matrix that contains all the background model for the inspected pixel, 
  % replacement for the the very slow x = squeeze( im(i,j,:) );
  Mymu=zeros(cfg.K,6);
  x = reshape( im(i,j,:), layers, 1 );
  i2=i-halfH;
  j2=j-halfW;
  oldj=round(r_h*j2*cos(rot_ang)-r_h*i2*sin(rot_ang)+halfW)+(round(dx));
  oldi=round(r_h*i2*cos(rot_ang)+r_h*j2*sin(rot_ang)+halfH)+(round(dy));
  
  if oldj<1 || oldi<1 || oldj>Wim || oldi>Him || mask(oldi, oldj)==1
       %siz = size( model(i,j,:,:) );
       for ks=1:cfg.K
        Mymu(ks,1:6)=[Mymodel(ks).mus(i,j,1) Mymodel(ks).mus(i,j,2) Mymodel(ks).mus(i,j,3)...
        Mymodel(ks).parameters(i,j,1) Mymodel(ks).parameters(i,j,2)  Mymodel(ks).parameters(i,j,3)];
       end  
       m=Mymu;
       %m = reshape( model(i,j,:,:), siz(3), siz(4) );
       m = newgauss2( x', m, idx, cfg ); 
       [idxmatch,no_match] = findmatch( x, m, idx, cfg,onemat );
  else
   
     % The following two lines replace the very slow m = squeeze( model(i,j,:,:) );
      %siz = size( model(oldi,oldj,:,:) );
      for ks=1:cfg.K
       Mymu(ks,1:6)=[Mymodel(ks).mus(oldi,oldj,1) Mymodel(ks).mus(oldi,oldj,2) Mymodel(ks).mus(oldi,oldj,3)...
          Mymodel(ks).parameters(oldi,oldj,1) Mymodel(ks).parameters(oldi,oldj,2)  Mymodel(ks).parameters(oldi,oldj,3)];
      end  
      m=Mymu;
     % m = reshape( model(oldi,oldj,:,:), siz(3), siz(4) );
      weights = m(:,idx.w);
     % Find the Gaussian that matches, see findmatch.
      [idxmatch,no_match] = findmatch( x, m, idx, cfg,onemat );

       
     % Initialize a new Gaussian if no match is found. The newgauss function
     % finds the Gaussian with minimal weight and replaces it by a new one.
      if no_match
        m = newgauss( x', m, idx, cfg );
     % If a match is found then decrease the weights of non-matching Gaussians
     % and re-normalize. Then update the model, see updatemodel.
      else
        % update weights
        % idx2update = find( [1:length(weights)]~=idxmatch );
        idx2update = ~idxmatch;
        weights(idx2update) = (1-cfg.alpha) * weights(idx2update);
        % renormalization
        weights = weights ./ sum(weights);
        m(:,idx.w) = weights;
        % update model
        m(idxmatch,:) = updatemodel( m(idxmatch,:), x', idx, cfg );
      end
  end
% After the update of the model we have to decide which Gaussians
% represent the background. It is important to note that the background
% may change due to lighting changes, object movement, etc. The set of 
% `background' Gaussians have to change accordingly.
  idxbck = findbckg( m, idx, cfg );
  m(:,idx.bckg) = 0;        % clear the labels
  m(idxbck,idx.bckg) = 1;   % label '1' means background
% Compose the background image. This part is actually not necessary for
% the segmentation. However, it is convenient to see how the background
% evolves over time: this composition can be moved outside the bckggm
% function. The background image is computed as a weighted average of the mean
% values of the Gaussian.
  if isscalar(idxbck)
    res.bckg(i,j,:) = m(idxbck,idx.mu);
  else % weighted average 
    w_rel = m(idxbck,idx.w); % related weights
    res.bckg(i,j,:) = sum(w_rel(:,ones(length(idx.mu),1)) .* ...
                      m(idxbck,idx.mu))/sum(m(idxbck,idx.w));
  end
    % res.bckg(i,j,:) = sum( repmat(m(idxbck,idx.w),1,length(idx.mu)) ...
    %                          .*m(idxbck,idx.mu) )./sum(m(idxbck,idx.w));
% The pixel is labeled as foreground if there was no match at all 
% or the matched Gaussian(s) does not belong to the background. 
  if sum(idxmatch)==0 | ~any(idxbck==find(idxmatch))
    res.segm(i,j) = 1;
  end
% Put the updated model matrix back in the data structure.
%Newmodel(i,j,:,:) = m;%%%**********************************
  for ks=1:cfg.K
      NewMymodel(ks).mus(i,j,:)=m(ks,1:3);
      NewMymodel(ks).parameters(i,j,:)=m(ks,4:6);
  end  
  end % for j
end % for i
%model=Newmodel;%%%****************************************
Mymodel=NewMymodel;
res.segm = logical(res.segm);

function [model,idx,Mymodel] = model_init2(im,cfg)
% Usage: [model,idx] = model_init(im,cfg)
% Function model_init initializes the Gaussians by using the 
% current observation and the predefined values stored in cfg.

r = size(im,1); % image height
c = size(im,2); % image width
b = size(im,3); % layers (3 for RGB, 2 for RG, 1 for intensity image)

% Initial parameters. The means are approximately uniformly
% distributed along the diagonal of the feature space. Initial
% variance is taken from cfg.var and the weights are set equally.
shift = cfg.thr * sqrt(cfg.var);
initmu = linspace( 0+shift, 1-shift, cfg.K-1 )';
initvars = cfg.var * ones(cfg.K-1,1);
initw = 1/cfg.K;

% Symbolic names for the indexes improve the readability of the code.
idx.mu = 1:b;       % means
idx.vars = b+1;     % variances
idx.w = idx.vars+1; % weights
idx.bckg = idx.w+1; % background flags

% Memory allocation for the complete model
model = zeros( r, c, cfg.K, idx.bckg );

Myinitmu=[0 initmu'];
Myinitbgr=[1 0 0 0];
for i=1:cfg.K
 Mymodel(i).mus=zeros(r,c,b);
 Mymodel(i).mus(:,:,:)=Myinitmu(i);
 Mymodel(i).parameters=zeros(r,c,3);
 Mymodel(i).parameters(:,:,1)=cfg.var;
 Mymodel(i).parameters(:,:,2)=initw;
 Mymodel(i).parameters(:,:,3)=Myinitbgr(i);
end
Mymodel(1).mus(:,:,:)=im;

return; % end of bckggm

function [model,idx,Mymodel] = model_init(im,cfg)
% Usage: [model,idx] = model_init(im,cfg)
% Function model_init initializes the Gaussians by using the 
% current observation and the predefined values stored in cfg.

r = size(im,1); % image height
c = size(im,2); % image width
b = size(im,3); % layers (3 for RGB, 2 for RG, 1 for intensity image)

% Initial parameters. The means are approximately uniformly
% distributed along the diagonal of the feature space. Initial
% variance is taken from cfg.var and the weights are set equally.
shift = cfg.thr * sqrt(cfg.var);
initmu = linspace( 0+shift, 1-shift, cfg.K-1 )';
initvars = cfg.var * ones(cfg.K-1,1);
initw = 1/cfg.K;

% Symbolic names for the indexes improve the readability of the code.
idx.mu = 1:b;       % means
idx.vars = b+1;     % variances
idx.w = idx.vars+1; % weights
idx.bckg = idx.w+1; % background flags

% Memory allocation for the complete model
model = zeros( r, c, cfg.K, idx.bckg );

Myinitmu=[0 initmu'];
Myinitbgr=[1 0 0 0];
for i=1:cfg.K
 Mymodel(i).mus=zeros(r,c,b);
 Mymodel(i).mus(:,:,:)=Myinitmu(i);
 Mymodel(i).parameters=zeros(r,c,3);
 Mymodel(i).parameters(:,:,1)=cfg.var;
 Mymodel(i).parameters(:,:,2)=initw;
 Mymodel(i).parameters(:,:,3)=Myinitbgr(i);
end
Mymodel(1).mus(:,:,:)=im;

% Matrix of initial parameters that is added to each pixel beside
% the current value
initpars = ...
  [initmu(:,ones(b,1)) initvars initw(ones(1,cfg.K-1),:) zeros(cfg.K-1,1)];

% For each pixel in the image take the current observation as initialization 
% of the  and assign initial values to the rest of the Gaussians
for i = 1:r
  for j = 1:c
    model(i,j,:,:) = [squeeze(im(i,j,:))' cfg.var initw 1; initpars];
  end
end
return; % end of bckggm

function [idxmatch,no_match] = findmatch(x,model,idx,cfg,onemat)
% Usage: [idxmatch,no_match] = findmatch(x,model,idx,cfg)
% findmatch compares the actual observation x
% with all Gaussians and selects the one that matches. If no match is
% found, an empty array is returned.
% 
% Absolute distances between the means and the current observation. 
%d = abs(model(:,idx.mu)' - x(:,ones(cfg.K,1)));
d = abs(model(:,[1 2 3])' - x(:,onemat)); % it is [1 2 3] because i have 3 colour channels

% Thresholding: Absolute distance in each channel (layer) is compared to
% cfg.thr-multiple of sigma-s. Only those Gaussians
% where distances in  layers are smaller than thresholds
% are accepted.
no_match = 0;
stds = sqrt( model(:,idx.vars) )';
thr_mat = cfg.thr * stds( ones(3,1), : );
id = find( all(d<thr_mat) );

% It can happen that several Gaussians match, and several possible strategies
% for selecting the `best' exist. Here one can select between two:
% either the one with minimal variance or the closest. Depending on time
% constraints, more advanced metrics such as Mahalanobis distance can be
% chosen.
if length(id)>1
  [minvar,idmin] = min( model(id,idx.vars) ); % minimal variance
  % [mind,idmin] = min( sum(d(:,id)) );       % closest 
  id = id(idmin);
elseif isempty(id)
  no_match = 1;
end
idxmatch(1:cfg.K) = false;
idxmatch(id) = 1;
return; % end of find match

function m = newgauss(x,m,idx,cfg)
% Usage: m = newgauss(x,m,idx,cfg)
% newgauss replaces the Gaussian with 
% the lowest weight by a new one. The current observation
% initializes the mean. The variance of the new Gaussian
% is computed as double that of the largest among the remaining
% Gaussians. The weight is the minimal one, and the weights 
% are re-normalized.
% 
[minw,idmin] = min(m(:,idx.w));           % minimal weight
m(idmin,idx.mu) = x;                      % assign new observation
m(idmin,idx.vars) = 2*max(m(:,idx.vars)); % twice the highest variance
m(idmin,idx.w) = minw;                    % take the minimal weight
m(:,idx.w) = m(:,idx.w)/sum(m(:,idx.w));  % re-normalize wights
return; % end of newgauss

function m = newgauss2(x,m,idx,cfg)
% Usage:when new pixel is introduces into the background
% 
[minw,idmin] = min(m(:,idx.w));           % minimal weight
[maxw,idmax] = max(m(:,idx.w)); 
m(idmin,idx.mu) = x;                      % assign new observation
m(idmin,idx.vars) = max(m(:,idx.vars)); % twice the highest variance*********************************************modified
m(idmin,idx.w) = maxw;                    % take the minimal weight
m(:,idx.w) = m(:,idx.w)/sum(m(:,idx.w));  % re-normalize wights
return; % end of newgauss


function m = updatemodel(m,x,idx,cfg)
% Usage: m = updatemodel(m,x,idx,cfg)
% First weight the learning constant by the `credibility' of the
% observation x. Clearly, the further the observation is from the 
% mean of the Gaussian the less influence it should have on the update.
% The probability is computed in a simplified way: the covariance
% matrix is assumed to be diagonal and all variances equal. It suffices
% to compute probability in one dimension. We skip the normalization
% of the probability value
% in order to be equal to 1 if x matches . 
% This speed-up learning trick is questionable.
% The probability value for the exact match is naturally lower for
% higher variance. This lower weighting of `flat' densities may be
% useful in some applications.
ro = cfg.alpha;% * exp( -0.5*(x(1)-m(idx.mu(1)))./sqrt(m(idx.vars)).^2 ); 
% update mean
m(idx.mu) = (1-ro)*m(idx.mu) + ro*x;
% update variances
%m(idx.vars) = (1-ro)*m(idx.vars) + ro*sum((m(:,idx.mu)-x).^2);
m(idx.vars) = max((1-ro)*m(idx.vars) + ro*sum((m(:,idx.mu)-x).^2),cfg.minvar);
return % end of updatemodel



function idxbck = findbckg(m,idx,cfg)
% Usage: idxbck = findbckg(m,idx,cfg)
% This selects the Gaussians that most probable represent background.
% The underlying idea is twofold: The background Gaussians presumably have 
% high weights which correspond to frequent observation, and
% low variances. The assumption is that we assume the total frequency
% of background is higher than cfg.T.
%
% First, sort Gaussians according to omega/sigma^2, best first.
sortcr = m(:,idx.w) ./ m(:,idx.vars);
[foo,idxsort] = sort( sortcr );
idxsort = idxsort(end:-1:1);
% Find the minimum number of sorted Gaussians whose sum
% exceeds the threshold cfg.T
idxbck = find( cumsum(m(idxsort,idx.w))>cfg.T );
B = idxbck(1);
idxbck = idxsort(1:B);
if isempty(idxbck)
  [foo,idxbck] = max( m(:,idx.w) );
end
return % end of findbckg

