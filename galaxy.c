#include <stdlib.h>
#include <math.h>

float rand_float_from_to(float, float);

typedef struct {
    double x, y, z;
} VECTOR;

typedef struct {
    VECTOR position;
    VECTOR acceleration;
    VECTOR velocity;
    double mass;
} STAR;

typedef struct {
    STAR stars[1000];
} GALAXY;

VECTOR star_position(double, double, double);

double Vector_magnitude(VECTOR);

GALAXY create_galaxy(float heightMagnitude, float heightFrequency, int numStar, VECTOR galaxy_center, VECTOR velocity) {

    GALAXY galaxy;

    for (int i = 0; i < numStar; i++) {
        //generovanie pozicie hviezdy, vrati vektor s hodnotami x,y,z
        VECTOR v = star_position(galaxy_center.x, galaxy_center.y, galaxy_center.z);

        v.y -= heightMagnitude * sin(Vector_magnitude(v) * heightFrequency);

        //nacitanie pozicie hviezdy
        double size = 1e7;
        galaxy.stars[i].position.x = v.x * size - size;
        galaxy.stars[i].position.y = v.y * size - size;
        galaxy.stars[i].position.z = v.z * size - size;

        //nastavenie pociatocnej rychlosti
        galaxy.stars[i].velocity.x = velocity.x;
        galaxy.stars[i].velocity.y = velocity.y;
        galaxy.stars[i].velocity.z = velocity.z;

        galaxy.stars[i].acceleration.x = 0;
        galaxy.stars[i].acceleration.y = 0;
        galaxy.stars[i].acceleration.z = 0;

        galaxy.stars[i].mass = 1e26;
    }
    return galaxy;
}


VECTOR star_position(double galaxy_center_x, double galaxy_center_y, double galaxy_center_z) {

    VECTOR v;
    int numArms = 2;
    double arm_separation_distance = 2 * 3.14 / numArms;
    float max_arm_offset = 1.0f;
    double rotationFactor = 10.0;

    // nastavenie vzdialenosti od stredu
    double distance = rand_float_from_to(0, 0.3f);
    distance = pow(distance, 1);

    // nastavenie uhla
    double angle = rand_float_from_to(0, 1) * 2 * M_PI;
    double arm_offset = rand_float_from_to(0, 0.3f) * max_arm_offset;

    arm_offset = arm_offset - max_arm_offset / 2;
    arm_offset = arm_offset * (1 / distance);

    double squared_arm_offset = pow(arm_offset, 2);
    if (arm_offset < 0)
        squared_arm_offset = squared_arm_offset * -1;
    arm_offset = squared_arm_offset;

    double rotation = distance * rotationFactor;

    angle = (int) (angle / arm_separation_distance) * arm_separation_distance + arm_offset + rotation;

    v.x = (cos(angle) * distance) + galaxy_center_x;
    v.z = (sin(angle) * distance) + galaxy_center_z;
    v.y = 0 + galaxy_center_y;

    return v;
}


float rand_float_from_to(float min, float max) {
    float result = (float) rand() / (float) (RAND_MAX) * (max - min) + min;
    return result;
}

double Vector_magnitude(VECTOR p) {
    double w = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
    return w;
}