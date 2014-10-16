module GreatCircle

export great_distance, great_circle

two_pi = 2.0 * pi


function great_circle(distance, azimuth, latitude, longitude, rmajor=6378137.0, rminor=6356752.3142)
    """
        Named arguments:
        distance = distance to travel, or numpy array of distances
        azimuth = angle, in DEGREES of HEADING from NORTH, or numpy array of azimuths
        latitude = latitude, in DECIMAL DEGREES, or numpy array of latitudes
        longitude = longitude, in DECIMAL DEGREES, or numpy array of longitudes
        rmajor = radius of earth's major axis. default=6378137.0 (WGS84)
        rminor = radius of earth's minor axis. default=6356752.3142 (WGS84)

        Returns a dictionary with:
        'latitude' in decimal degrees
        'longitude' in decimal degrees
        'reverse_azimuth' in decimal degrees
    """

    azimuth = deg2rad(azimuth)
    latitude = deg2rad(latitude)
    longitude = deg2rad(longitude)
    f = (rmajor - rminor) / rmajor

    new_lat, new_lon, reverse_azimuth = vincentypt(f, rmajor, latitude, longitude, azimuth, distance)
    {"latitude"=>rad2deg(new_lat), "longitude"=>rad2deg(new_lon), "reverse_azimuth"=>rad2deg(reverse_azimuth)}
end


function great_distance(start_latitude, start_longitude, end_latitude, end_longitude, rmajor=6378137.0, rminor=6356752.3142)
    """
        Named arguments:
        start_latitude = starting latitude, in DECIMAL DEGREES
        start_longitude = starting longitude, in DECIMAL DEGREES
        end_latitude = ending latitude, in DECIMAL DEGREES
        end_longitude = ending longitude, in D|ECIMAL DEGREES
        rmajor = radius of earth's major axis. default=6378137.0 (WGS84)
        rminor = radius of earth's minor axis. default=6356752.3142 (WGS84)

        Returns a dictionaty with:
        'distance' in meters
        'azimuth' in decimal degrees
        'reverse_azimuth' in decimal degrees
    """

    start_latitude = deg2rad(start_latitude)
    start_longitude = deg2rad(start_longitude)
    end_latitude = deg2rad(end_latitude)
    end_longitude = deg2rad(end_longitude)
    f = (rmajor - rminor) / rmajor

    distance, angle, reverse_angle = vincentydist(f, rmajor, start_latitude, start_longitude, end_latitude, end_longitude)
    {"distance"=>rad2deg(distance), "azimuth"=>rad2deg(angle), "reverse_azimuth"=>rad2deg(reverse_angle)}
end


# -----------------------------------------------------------------------
# | Algrothims from Geocentric Datum of Australia Technical Manual |
# | |
# | http://www.anzlic.org.au/icsm/gdatum/chapter4.html |
# | |
# | This page last updated 11 May 1999 |
# | |
# | Computations on the Ellipsoid |
# | |
# | There are a number of formulae that are available |
# | to calculate accurate geodetic positions, |
# | azimuths and distances on the ellipsoid. |
# | |
# | Vincenty's formulae (Vincenty, 1975) may be used |
# | for lines ranging from a few cm to nearly 20,000 km, |
# | with millimetre accuracy. |
# | The formulae have been extensively tested |
# | for the Australian region, by comparison with results |
# | from other formulae (Rainsford, 1955 & Sodano, 1965). |
# | |
# | * Inverse problem: azimuth and distance from known |
# | latitudes and longitudes |
# | * Direct problem: Latitude and longitude from known |
# | position, azimuth and distance. |
# | * Sample data |
# | * Excel spreadsheet |
# | |
# | Vincenty's Inverse formulae |
# | Given: latitude and longitude of two points |
# | (phi1, lembda1 and phi2, lembda2), |
# | Calculate: the ellipsoidal distance (s) and |
# | forward and reverse azimuths between the points (alpha12, alpha21). |
# | |
# -----------------------------------------------------------------------
function vincentypt(f::Number, a::Number, phi1::Number, lembda1::Number, alpha12::Number, s::Number)
    """
    Returns: lat and long of projected point and reverse azimuth,
    given a reference point and a distance and azimuth to project.
    lats, longs and azimuths are passed in RADIANS

    Returns ( phi2, lambda2, alpha21 ) as a tuple, all in radians
    """
    #----------------------------------------------------------------------------
    # Vincenty's Direct formulae |
    # Given: latitude and longitude of a point (phi1, lembda1) and |
    # the geodetic azimuth (alpha12) |
    # and ellipsoidal distance in metres (s) to a second point, |
    # |
    # Calculate: the latitude and longitude of the second point (phi2, lembda2) |
    # and the reverse azimuth (alpha21). |
    # |
    #----------------------------------------------------------------------------
    if alpha12 < 0.0
        alpha12 = alpha12 + two_pi
    end
    if alpha12 > two_pi
        alpha12 = alpha12 - two_pi
    end

    b = a .* (1.0 - f)

    TanU1 = (1 - f) .* tan(phi1)
    U1 = atan( TanU1 )
    sigma1 = atan2( TanU1, cos(alpha12) )
    Sinalpha = cos(U1) .* sin(alpha12)
    cosalpha_sq = 1.0 - Sinalpha .* Sinalpha

    u2 = cosalpha_sq .* (a .* a - b .* b ) ./ (b .* b)
    A = 1.0 + (u2 ./ 16384) .* (4096 + u2 .* (-768 + u2 .* (320 - 175 .* u2) ) )
    B = (u2 ./ 1024) .* (256 + u2 .* (-128 + u2 .* (74 - 47 .* u2) ) )

    # Starting with the approximation
    sigma = (s ./ (b .* A))

    # Not moving anywhere. We can return the location that was passed in.
    if sigma == 0
        return phi1, lembda1, alpha12
    end

    last_sigma = 2.0 .* sigma + 2.0 # something impossible

    # Iterate the following three equations
    # until there is no significant change in sigma
    # two_sigma_m , delta_sigma
    two_sigma_m = 0.
    while ( abs( (last_sigma - sigma) ./ sigma) > 1.0e-9 )
        two_sigma_m = 2 .* sigma1 + sigma
        delta_sigma = B .* sin(sigma) .* ( cos(two_sigma_m) + (B./4) .* (cos(sigma) .* (-1 + 2 .* cos(two_sigma_m).^2 - (B./6) .* cos(two_sigma_m) .* (-3 + 4 .* sin(sigma).^2) .* (-3 + 4 .* cos(two_sigma_m).^2 ))))
        last_sigma = sigma
        sigma = (s ./ (b .* A)) + delta_sigma
    end

    phi2 = atan2 ( (sin(U1) .* cos(sigma) + cos(U1) .* sin(sigma) .* cos(alpha12) ), ((1-f) .* sqrt( Sinalpha.^2 + (sin(U1) .* sin(sigma) - cos(U1) .* cos(sigma) .* cos(alpha12)).^2)))

    lembda = atan2( (sin(sigma) .* sin(alpha12 )), (cos(U1) .* cos(sigma) - sin(U1) .* sin(sigma) .* cos(alpha12)))

    C = (f./16) .* cosalpha_sq .* (4 + f .* (4 - 3 .* cosalpha_sq ))

    omega = lembda - (1-C) .* f .* Sinalpha .* (sigma + C .* sin(sigma) .* (cos(two_sigma_m) + C .* cos(sigma) .* (-1 + 2 .* cos(two_sigma_m).^2 )))

    lembda2 = lembda1 + omega

    alpha21 = atan2 ( Sinalpha, (-sin(U1) .* sin(sigma) + cos(U1) .* cos(sigma) .* cos(alpha12)))
    alpha21 = alpha21 + two_pi ./ 2.0

    if alpha21 < 0.0
        alpha21 = alpha21 + two_pi
    end
    if alpha21 > two_pi
        alpha21 = alpha21 - two_pi
    end
    return phi2, lembda2, alpha21

end

function vincentypt(f::Number, a::Number, phi1::Array, lembda1::Array, alpha12::Array, s::Array)
    lenphi = length(phi1)
    @assert lenphi == length(lembda1)
    @assert lenphi == length(alpha12)
    @assert lenphi == length(s)
    phi2 = zeros(phi1)
    lembda2 = zeros(phi1)
    alpha21 = zeros(phi1)
    for i=1:lenphi
        @inbounds phi2[i], lembda2[i], alpha21[i] = vincentypt(f, a, phi1[i], lembda1[i], alpha12[i], s[i])
    end
    return phi2, lembda2, alpha21
end
# f, rmajor, latitude, longitude, azimuth, distance
vincentypt(f::Number, a::Number, phi1::Array, lembda1::Array, alpha12::Array, s::Number) = vincentypt(f, a, phi1, lembda1, alpha12, s*ones(lembda1))
vincentypt(f::Number, a::Number, phi1::Array, lembda1::Array, alpha12::Number, s::Array) = vincentypt(f, a, phi1, lembda1, alpha12*ones(lembda1), s)
vincentypt(f::Number, a::Number, phi1::Array, lembda1::Array, alpha12::Number, s::Number) = vincentypt(f, a, phi1, lembda1, alpha12*ones(lembda1), s*ones(lembda1))


function vincentydist(f::Number, a::Number, phi1::Number, lembda1::Number, phi2::Number, lembda2::Number)
    """
    Returns the distance between two geographic points on the ellipsoid
    and the forward and reverse azimuths between these points.
    lats, longs and azimuths are in radians, distance in meters

    Returns ( s, alpha12, alpha21 ) as a tuple
    """

    if (abs( phi2 - phi1 ) < 1e-8) & ( abs( lembda2 - lembda1) < 1e-8 )
        return 0.0, 0.0, 0.0
    end

    b = a .* (1.0 - f)

    TanU1 = (1 - f) .* tan( phi1 )
    TanU2 = (1 - f) .* tan( phi2 )

    U1 = atan(TanU1)
    U2 = atan(TanU2)

    lembda = lembda2 - lembda1
    last_lembda = -4000000.0 # an impossibe value
    omega = lembda

    # Iterate the following equations,
    # until there is no significant change in lembda
    while ( last_lembda < -3000000.0 | lembda != 0) & (abs( (last_lembda - lembda)./lembda) > 1.0e-9 )
        sqr_sin_sigma = (cos(U2) .* sin(lembda)).^2 + ( (cos(U1) .* sin(U2) - sin(U1) .* cos(U2) .* cos(lembda) )).^2
        Sin_sigma = sqrt( sqr_sin_sigma )
        Cos_sigma = sin(U1) .* sin(U2) + cos(U1) .* cos(U2) .* cos(lembda)
        sigma = atan2( Sin_sigma, Cos_sigma )
        Sin_alpha = cos(U1) .* cos(U2) .* sin(lembda) ./ sin(sigma)
        alpha = asin( Sin_alpha )
        Cos2sigma_m = cos(sigma) - (2 .* sin(U1) .* sin(U2) ./ cos(alpha).^2 )
        C = (f./16) .* cos(alpha).^2 .* (4 + f .* (4 - 3 .* cos(alpha).^2))
        last_lembda = lembda
        lembda = omega + (1-C) .* f .* sin(alpha) .* (sigma + C .* sin(sigma) .* (Cos2sigma_m + C .* cos(sigma) .* (-1 + 2 .* Cos2sigma_m.^2 )))
    end

    u2 = cos(alpha).^2 .* (a.*a-b.*b) ./ (b.*b)

    A = 1 + (u2./16384) .* (4096 + u2 .* (-768 + u2 .* (320 - 175 .* u2)))

    B = (u2./1024) .* (256 + u2 .* (-128+ u2 .* (74 - 47 .* u2)))

    delta_sigma = B .* Sin_sigma .* (Cos2sigma_m + (B./4) .* (Cos_sigma .* (-1 + 2 .* Cos2sigma_m.^2 ) - (B./6) .* Cos2sigma_m .* (-3 + 4 .* sqr_sin_sigma) .* (-3 + 4 .* Cos2sigma_m.^2 )))

    s = b .* A .* (sigma - delta_sigma)

    alpha12 = atan2( (cos(U2) .* sin(lembda)), (cos(U1) .* sin(U2) - sin(U1) .* cos(U2) .* cos(lembda)))
    alpha21 = atan2( (cos(U1) .* sin(lembda)), (-sin(U1) .* cos(U2) + cos(U1) .* sin(U2) .* cos(lembda)))

    if alpha12 < 0.0
        alpha12 = alpha12 + two_pi
    end
    if alpha12 > two_pi
        alpha12 = alpha12 - two_pi
    end
    alpha21 = alpha21 + two_pi ./ 2.0

    if alpha21 < 0.0
        alpha21 = alpha21 + two_pi
    end
    if alpha21 > two_pi
        alpha21 = alpha21 - two_pi
    end
    return s, alpha12, alpha21
end

function vincentydist(f::Number, a::Number, phi1::Array, lembda1::Array, phi2::Array, lembda2::Array)
    lenphi = length(phi1)
    @assert lenphi == length(lembda1)
    @assert lenphi == length(phi2)
    @assert lenphi == length(lembda2)
    s = zeros(phi1)
    alpha12 = zeros(phi1)
    alpha21 = zeros(phi1)
    for i=1:lenphi
        @inbounds s[i], alpha12[i], alpha21[i] = vincentydist(f, a, phi1[i], lembda1[i], phi2[i], lembda2[i])
    end
    return s, alpha12, alpha21
end


end # module
