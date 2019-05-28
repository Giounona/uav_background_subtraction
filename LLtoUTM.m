%this function converts the lattitude and longitude to easting and northing
%values

%input: lat=lattitude
%       lon=longitude
%        
%output: east=easting
%        north=northing
%        zoneNumber 

%this source file is the copyright property of the Kingston University (see README.txt file)
%Author: Giounona Tzanidou
%date: 3-07-2015



function [east north zoneNumber ]=LLtoUTM(lat, lon)
UTMScaleFactor = 0.9996;
C_WGS84_a = 6378137.000;
C_WGS84_b = 6356752.31414;
M_PI=pi;

	zone = round(((lon + 180.0)/6.0) + 1.0);
    
	if ( (lat >= 56.0) && (lat < 64.0) && (lon >= 3.0) && (lon < 12.0) )
		zone = 32;
    end
    
	if ( (lat >= 72.0) && (lat < 84.0) )
	
		if		( (lon >= 0.0) && (lon < 9.0) )
                zone = 31;
        elseif ( (lon >= 9.0) && (lon < 21.0) ) 
                zone = 33;
        elseif ( (lon >= 21.0) && (lon < 33.0) )
                zone = 35;
        elseif ( (lon >= 33.0) && (lon < 42.0) ) 
                zone = 37;    
        end
    end

	phi = lat / 180 * M_PI;
	lambda = lon / 180 * M_PI;
	lambda0 = (-183.0 + (zone * 6.0)) / 180.0 * M_PI;

	%/* Precalculate ep2 */
	ep2 = ((C_WGS84_a^2.0) - (C_WGS84_b^ 2.0)) / (C_WGS84_b^ 2.0);

	%/* Precalculate nu2 */
	nu2 = ep2 * (cos(phi)^ 2.0);

	%/* Precalculate N */
	N = (C_WGS84_a^ 2.0) / (C_WGS84_b * sqrt(1.0 + nu2));

	%/* Precalculate t */
	t = tan(phi);
	t2 = t * t;
	tmp = (t2 * t2 * t2) - (t^ 6.0);

	%/* Precalculate l */
	l = lambda - lambda0;

	%/* Precalculate coefficients for l**n in the equations below
	%so a normal human being can read the expressions for easting
	%and northing
	%-- l**1 and l**2 have coefficients of 1.0 */
    
	l3coef = 1.0 - t2 + nu2;
	l4coef = 5.0 - t2 + 9 * nu2 + 4.0 * (nu2 * nu2);
	l5coef = 5.0 - 18.0 * t2 + (t2 * t2) + 14.0 * nu2 - 58.0 * t2 * nu2;
	l6coef = 61.0 - 58.0 * t2 + (t2 * t2) + 270.0 * nu2 - 330.0 * t2 * nu2;
	l7coef = 61.0 - 479.0 * t2 + 179.0 * (t2 * t2) - (t2 * t2 * t2);
    l8coef = 1385.0 - 3111.0 * t2 + 543.0 * (t2 * t2) - (t2 * t2 * t2);

	%/* Calculate easting (x) */
	 xy0 = N * cos(phi) * l+ (N / 6.0 * (cos(phi)^ 3.0) * l3coef * (l^ 3.0))+ (N / 120.0 * (cos(phi)^ 5.0) * l5coef * (l^ 5.0))+ (N / 5040.0 * (cos(phi)^ 7.0) * l7coef * (l^ 7.0));

	%/* Calculate northing (y) */

	%/* Precalculate n */
	 n = (C_WGS84_a - C_WGS84_b) / (C_WGS84_a + C_WGS84_b);

	%/* Precalculate alpha */
	 alpha = ((C_WGS84_a + C_WGS84_b) / 2.0) * (1.0 + ((n^ 2.0) / 4.0) + ((n^ 4.0) / 64.0));

	%/* Precalculate beta */
	 beta = (-3.0 * n / 2.0) + (9.0 * (n^ 3.0) / 16.0) + (-3.0 * (n^ 5.0) / 32.0);

	%/* Precalculate gamma */
	 gamma = (15.0 * (n^ 2.0) / 16.0) + (-15.0 * (n^ 4.0) / 32.0);

	%/* Precalculate delta */
	 delta = (-35.0 * (n^ 3.0) / 48.0) + (105.0 * (n^ 5.0) / 256.0);

	%/* Precalculate epsilon */
	 epsilon = (315.0 * (n^ 4.0) / 512.0);

	%/* Now calculate the sum of the series and return */
	 ArcLengthOfMeridian = alpha* (phi + (beta * sin(2.0 * phi))+ (gamma * sin(4.0 * phi))+ (delta * sin(6.0 * phi))+ (epsilon * sin(8.0 * phi)));

	 xy1 = ArcLengthOfMeridian+ (t / 2.0 * N * (cos(phi)^ 2.0) * (l^ 2.0))+ (t / 24.0 * N * (cos(phi)^ 4.0) * l4coef * (l^ 4.0))+ (t / 720.0 * N * (cos(phi)^ 6.0) * l6coef * (l^ 6.0))+ (t / 40320.0 * N * (cos(phi)^ 8.0) * l8coef * (l^ 8.0));

	%/* Adjust easting and northing for UTM system. */
	xy0 = xy0 * UTMScaleFactor + 500000.0;
	xy1 = xy1 * UTMScaleFactor;
	if xy1 < 0.0
        xy1 = xy1 + 10000000.0;
    end

	east = xy0;
	north = xy1;
	zoneNumber = zone;
	if (lat >= 0.0) 
         hemisphere = 'N'; 
    else
         hemisphere = 'S'; 
    end
            
    
