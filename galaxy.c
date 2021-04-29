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
    STAR stars[2000];
    VECTOR center;
    double mass;
} GALAXY;

VECTOR star_position(double, double, double);

double Vector_magnitude(VECTOR);

VECTOR norm_vector(VECTOR);

VECTOR cross_vector(VECTOR, VECTOR);

double orbital_vel(double, double);

GALAXY create_galaxy(float heightMagnitude, float heightFrequency, int numStar, VECTOR galaxy_center, VECTOR velocity) {

    GALAXY galaxy;

    galaxy.center.x = galaxy_center.x;
    galaxy.center.y = galaxy_center.y;
    galaxy.center.z = galaxy_center.z;
    for (int i = 0; i < numStar; i++) {
        //generovanie pozicie hviezdy, vrati vektor s hodnotami x,y,z
        VECTOR v = star_position(galaxy_center.x, galaxy_center.y, galaxy_center.z);

        v.y -= heightMagnitude * sin(Vector_magnitude(v) * heightFrequency);

        //nacitanie pozicie hviezdy
        double size = 1e7;
        galaxy.stars[i].position.x = (v.x * size - size) ;
        galaxy.stars[i].position.y = (v.y * size - size) ;
        galaxy.stars[i].position.z = (v.z * size - size) ;
        //pusunutie na center

        galaxy.stars[i].position.x +=size;
        galaxy.stars[i].position.y +=size;
        galaxy.stars[i].position.z +=size;
        galaxy.stars[i].mass = 1e24;
        galaxy.mass = galaxy.stars[i].mass;
        //nastavenie pociatocnej rychlosti
        VECTOR up;
        up.x=0;
        up.y=1;
        up.z=0;
        VECTOR vec = cross_vector(galaxy.stars[i].position,up);
        VECTOR vec1 = norm_vector(vec);
        VECTOR relative_vel;
        double orbital_velocity;
        //vzdialenost od stredu
        VECTOR vz;
        vz.x = galaxy.center.x - galaxy.stars[i].position.x;
        vz.y = galaxy.center.y - galaxy.stars[i].position.y;
        vz.z = galaxy.center.z - galaxy.stars[i].position.z;
        orbital_velocity = orbital_vel(galaxy.stars[i].mass,Vector_magnitude(galaxy.stars[i].position));
        relative_vel.x = (vec1.x * orbital_velocity)*0.5;
        relative_vel.y = (vec1.y * orbital_velocity)*0.5;
        relative_vel.z = (vec1.z * orbital_velocity)*0.5;


        galaxy.stars[i].velocity.x = velocity.x + relative_vel.x;
        galaxy.stars[i].velocity.y = velocity.y + relative_vel.y;
        galaxy.stars[i].velocity.z = velocity.z + relative_vel.z;

        galaxy.stars[i].acceleration.x = 0;
        galaxy.stars[i].acceleration.y = 0;
        galaxy.stars[i].acceleration.z = 0;


    }

    galaxy.mass = galaxy.mass * 100;
    return galaxy;
}


VECTOR star_position(double galaxy_center_x, double galaxy_center_y, double galaxy_center_z) {

    VECTOR v;
    int numArms = 2;
    double arm_separation_distance = 2 * 3.14 / numArms;
    float max_arm_offset = 1.5f;
    double rotationFactor = 7.0;

    // nastavenie vzdialenosti od stredu
    double distance = rand_float_from_to(0, 0.7f);
    distance = pow(distance, 1);

    // nastavenie uhla
    double angle = rand_float_from_to(0, 1) * 2 * M_PI;
    double arm_offset = rand_float_from_to(0, 0.7f) * max_arm_offset;

    arm_offset = arm_offset - max_arm_offset / 2;
    arm_offset = arm_offset * (1 / distance);

    double squared_arm_offset = pow(arm_offset, 2);
    if (arm_offset < 0)
        squared_arm_offset = squared_arm_offset * -1;
    arm_offset = squared_arm_offset;

    double rotation = distance * rotationFactor;

    angle = (int) (angle / arm_separation_distance) * arm_separation_distance + arm_offset + rotation;

    v.x = (cos(angle) * distance) ;
    v.z = (sin(angle) * distance) ;
    v.y = 0 ;

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

VECTOR norm_vector(VECTOR vector){
    VECTOR vector1;
    double magnitude = Vector_magnitude(vector);
    double smag = sqrt(magnitude);
    if(magnitude>0){
        vector1.x = vector.x/smag;
        vector1.y = vector.y/smag;
        vector1.z = vector.z/smag;
    } else {
        vector1.x = 0;
        vector1.y = 0;
        vector1.z = 0;
    }
    return vector1;
}

VECTOR cross_vector(VECTOR first, VECTOR second){
    VECTOR out_vector;
    out_vector.x = first.y * second.z - first.z * second.y;
    out_vector.y = first.z * second.x - first.x * second.z;
    out_vector.z = first.x * second.y - first.y * second.x;

    return out_vector;
}

double orbital_vel(double mass,double radius){
    double G = 6.6742367e-11; // m^3.kg^-1.s^-2
    return sqrt(G*mass/radius);
}
