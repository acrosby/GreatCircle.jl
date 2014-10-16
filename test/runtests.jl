using GreatCircle
using Base.Test

# One decimal degree is 111000m
latitude = 40.0
longitude = -76.0

azimuth = 90
new_gc = great_circle(111000, azimuth, latitude, longitude)
# We should have gone to the right
@test new_gc["longitude"] > longitude + 0.9

azimuth = 270
new_gc = great_circle(111000, azimuth, latitude, longitude)
# We should have gone to the left
@test new_gc["longitude"] < longitude - 0.9

azimuth = 180
new_gc = great_circle(111000, azimuth, latitude, longitude)
# We should have gone down
@test new_gc["latitude"] < latitude - 0.9

azimuth = 0
new_gc = great_circle(111000, azimuth, latitude, longitude)
# We should have gone up
@test new_gc["latitude"] > latitude + 0.9

azimuth = 315
new_gc = great_circle(111000, azimuth, latitude, longitude)
# We should have gone up and to the left
@test new_gc["latitude"] > latitude + 0.45
@test new_gc["longitude"] < longitude - 0.45


# One decimal degree is 111000m
latitude = [40.0:5.:60.0]
longitude = [-40.0:-5.:-60.0]

azimuth = ones(longitude) * 90
new_gc = great_circle(111000, azimuth, latitude, longitude)
# We should have gone to the right
@test all(new_gc["longitude"] .> longitude + 0.9)

azimuth = 270
distance = 111000 * ones(longitude)
new_gc = great_circle(distance, azimuth, latitude, longitude)
# We should have gone to the left
@test all(new_gc["longitude"] .< longitude - 0.9)

azimuth = 180
new_gc = great_circle(111000, azimuth, latitude, longitude)
# We should have gone down
@test all(new_gc["latitude"] .< latitude - 0.9)

azimuth = 0
new_gc = great_circle(111000, azimuth, latitude, longitude)
# We should have gone up
@test all(new_gc["latitude"] .> latitude + 0.9)

azimuth = 315
new_gc = great_circle(111000, azimuth, latitude, longitude)
# We should have gone up and to the left
@test all(new_gc["latitude"] .> latitude + 0.45)
@test all(new_gc["longitude"] .< longitude - 0.45)
