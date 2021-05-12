#include <GL/glut.h>
#include <GLFW/glfw3.h>
#include "galaxy.c"
#include "time.h"
#include "stdio.h"
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>

void reshape(int, int);

void myInit();

void myDraw();

void gravity_calculate_acceleration(int,int);

void gravity_calculate_acceleration2(int,int);

void *calculate_galaxy2(void*);

void* calculate_galaxy1(void*);

void start_cal();

void cal_velocity(int);

#define num_star 2000

double half_time_step;
double time_step;

GALAXY galaxy, galaxy2;

int main(int argc, char *argv[]) {

    srand((unsigned int) time(NULL));
    float height_frequency = rand_float_from_to(0, 1.5f);
    float height_magnitude = rand_float_from_to(0, 0.7f);

    //prva galaxia
    VECTOR galaxy_center, velocity;
    galaxy_center.x = 0.5;
    galaxy_center.y = 0;
    galaxy_center.z = 0;

   //velocity.x = 15 * 1e5;
    velocity.x = 0;
    velocity.y = 0;
    velocity.z = 0;


    galaxy = create_galaxy(height_magnitude, height_frequency, num_star, galaxy_center, velocity);

    //druha galaxia
//    VECTOR galaxy_center_2, velocity2;
//    galaxy_center_2.x = 1;
//    galaxy_center_2.y = 1;
//    galaxy_center_2.z = 1;
//
//    velocity2.x = -15 * 1e5;
//    velocity2.y = 0;
//    velocity2.z = 0;
//
//    galaxy2 = create_galaxy(height_magnitude, height_frequency, num_star, galaxy_center_2, velocity2);

    time_step = 1.0/64.0;
    half_time_step = 0.5 * time_step;

    /* Initialize window system */
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowSize(1920, 1080);
    glutCreateWindow("Galaxy");

    /* Initialize graphics */
    myInit();

    /* Display callback and enter event loop */
    glutDisplayFunc(myDraw);
    glutReshapeFunc(reshape);
    glutIdleFunc(myDraw);

    glutMainLoop();

    return 1;
}

void myInit() {
    /* Background color */
    glClearColor(0.0, 0.0, 0.0, 1.0);

    /* Projection */
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

}

void reshape(int w, int h) {
    if (h == 0)
        h = 1;
    double ratio = w * 1.0 / h;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, w, h);
    gluPerspective(45.0f, ratio, 0.1f * 1e7, 100.0f * 1e7);
    glMatrixMode(GL_MODELVIEW);
}

void myDraw() {

    /* Clear the screen */
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();

    gluLookAt(0.0f * 1e7, 3.0f * 1e7, -0.0f * 1e7, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
    glPointSize(1.0);
    glBegin(GL_POINTS);

    start_cal();

    glColor3f(1.0, 1.0, 1.0);
    for (int i = 0; i < num_star; i++) {
        glVertex3f(galaxy.stars[i].position.x, galaxy.stars[i].position.y, galaxy.stars[i].position.z);
    }
    glColor3f(1.0, 0.2, 0.2);
//    for (int i = 0; i < num_star; i++) {
//        glVertex3f(galaxy2.stars[i].position.x, galaxy2.stars[i].position.y, galaxy2.stars[i].position.z);
//    }
    glVertex3f(0, 0, 0);
    //sleep(1);
    glEnd();
    glutSwapBuffers();
}

void start_cal(){
    int range1[2], range2[2];
    void *send_range1, *send_range2;
    double *receive_range1, *receive_range2, *receive_range3, *receive_range4;
    range1[0]=0;
    range1[1]= num_star/2;
    range2[0]=num_star/2;
    range2[1]=num_star;

    send_range1 = &range1;
    send_range2 = &range2;

    pthread_t tid1, tid2, tid3, tid4;

    //prva galaxia
    pthread_create(&tid1, NULL, calculate_galaxy1, (void *) send_range1);//start thread
    pthread_create(&tid2, NULL, calculate_galaxy1, (void *) send_range2);
    //druha
//    pthread_create(&tid3, NULL, calculate_galaxy2, (void *) send_range1);
//    pthread_create(&tid4, NULL, calculate_galaxy2, (void *) send_range2);

    pthread_join(tid1, (void **) &receive_range1);
    pthread_join(tid2, (void **) &receive_range2);
//    pthread_join(tid3, (void **) &receive_range3);
//    pthread_join(tid4, (void **) &receive_range4);


}

void* calculate_galaxy1(void* param){
    int *ran = (int *) param;
    int start, end;

    start = *(ran);
    end = *(ran + 1);

    for(int i = start; i<end; i++){
        galaxy.stars[i].position.x += half_time_step * galaxy.stars[i].velocity.x;
        galaxy.stars[i].position.y += half_time_step * galaxy.stars[i].velocity.y;
        galaxy.stars[i].position.z += half_time_step * galaxy.stars[i].velocity.z;
    }

   // gravity_calculate_acceleration(start,end);

    for(int i = start; i<end; i++){

        galaxy.stars[i].velocity.x += time_step * galaxy.stars[i].acceleration.x;
        galaxy.stars[i].velocity.y += time_step * galaxy.stars[i].acceleration.y;
        galaxy.stars[i].velocity.z += time_step * galaxy.stars[i].acceleration.z;

     //   cal_velocity(i);
        galaxy.stars[i].position.x += half_time_step * galaxy.stars[i].velocity.x;
        galaxy.stars[i].position.y += half_time_step * galaxy.stars[i].velocity.y;
        galaxy.stars[i].position.z += half_time_step * galaxy.stars[i].velocity.z;
    }



}

void *calculate_galaxy2(void* param) {
    int *ran = (int *) param;
    int start, end;

    start = *(ran);
    end = *(ran + 1);

    for (int i = start; i < end; i++) {
        galaxy2.stars[i].position.x += half_time_step * galaxy2.stars[i].velocity.x;
        galaxy2.stars[i].position.y += half_time_step * galaxy2.stars[i].velocity.y;
        galaxy2.stars[i].position.z += half_time_step * galaxy2.stars[i].velocity.z;
    }


  //  gravity_calculate_acceleration2(start,end);
    for (int i = start; i < end; i++) {
        galaxy2.stars[i].velocity.x += time_step * galaxy2.stars[i].acceleration.x;
        galaxy2.stars[i].velocity.y += time_step * galaxy2.stars[i].acceleration.y;
        galaxy2.stars[i].velocity.z += time_step * galaxy2.stars[i].acceleration.z;
        cal_velocity(i);
        galaxy2.stars[i].position.x += half_time_step * galaxy2.stars[i].velocity.x;
        galaxy2.stars[i].position.y += half_time_step * galaxy2.stars[i].velocity.y;
        galaxy2.stars[i].position.z += half_time_step * galaxy2.stars[i].velocity.z;
    }
}

//vypocet pre prvu
void gravity_calculate_acceleration(int start, int end) {
    double G = 6.6742367e-11; // m^3.kg^-1.s^-2
    double EPS =4e5;
    double orbital_vel;
    double distance;
    VECTOR v,norm_v;
  //  double EPS =0;
    for (int i = start; i < end; i++) {
        galaxy.stars[i].acceleration.x = 0;
        galaxy.stars[i].acceleration.y = 0;
        galaxy.stars[i].acceleration.z = 0;
        for (int j = 0; j < num_star; j++) {
            if (j == i) {
                continue;
            }
//            vypocet vzdielonosti medzi hviezdami vramci svojej galaxie
            double dx = galaxy.stars[i].position.x - galaxy.stars[j].position.x;
            double dy = galaxy.stars[i].position.y - galaxy.stars[j].position.y;
            double dz = galaxy.stars[i].position.z - galaxy.stars[j].position.z;

            double dist = sqrt(dx * dx + dy * dy+ dz * dz);
            double preff = pow(dist,2) + pow(EPS,2);
            double pref = -G/pow(preff,1.5)*galaxy.stars[j].mass;
//            double pref = -G*galaxy.stars[j].mass/pow(dist,3);
            galaxy.stars[i].acceleration.x += pref * dx * 1e1;
            galaxy.stars[i].acceleration.y += pref * dy * 1e1;
            galaxy.stars[i].acceleration.z += pref * dz * 1e1;


//            double dxs = galaxy.center.x - galaxy.stars[i].position.x;
//            double dys = galaxy.center.y - galaxy.stars[i].position.y;
//            double dzs = galaxy.center.z - galaxy.stars[i].position.z;
//
//            double dists = sqrt(dxs * dxs + dys * dys+ dzs * dzs);
//            double preffs = pow(dists,2) + pow(EPS,2);
//            double prefs = G/pow(preffs,1.5)*galaxy.mass;
//            double prefs = G*galaxy.mass/pow(dists,3);
//            galaxy.stars[i].acceleration.x += prefs * dxs * 1e1;
//            galaxy.stars[i].acceleration.y += prefs * dys * 1e1;
//            galaxy.stars[i].acceleration.z += prefs * dzs * 1e1;

//            v.x = galaxy.stars[i].position.x - galaxy.stars[j].position.x;
//            v.y = galaxy.stars[i].position.y - galaxy.stars[j].position.y;
//            v.z = galaxy.stars[i].position.z - galaxy.stars[j].position.z;
//
//            double dis = v.x*v.x+v.y*v.y+v.z*v.z;
//            distance = sqrt(dis);
//            orbital_vel = -G*galaxy.stars[j].mass/distance/distance;
//            norm_v = norm_vector(v);
//            galaxy.stars[i].velocity.x += (norm_v.x*orbital_vel)*time_step;
//            galaxy.stars[i].velocity.y += (norm_v.y*orbital_vel)*time_step;
//            galaxy.stars[i].velocity.z += (norm_v.z*orbital_vel)*time_step;

        }
    }

}

//vypocet pre druhu
void gravity_calculate_acceleration2(int start, int end) {
    double G = 6.6742367e-11; // m^3.kg^-1.s^-2
    double EPS =3e4;
    for (int i = start; i < end; i++) {
        galaxy2.stars[i].acceleration.x = 0;
        galaxy2.stars[i].acceleration.y = 0;
        galaxy2.stars[i].acceleration.z = 0;
        double sum_x = 0, sum_y = 0, sum_z = 0;
        for (int j = 0; j < num_star; j++) {
            if (j == i) {
                double dx2 = galaxy2.stars[i].position.x - galaxy.stars[j].position.x;
                double dy2 = galaxy2.stars[i].position.y - galaxy.stars[j].position.y;
                double dz2 = galaxy2.stars[i].position.z - galaxy.stars[j].position.z;
                double dist2 = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
                double preff2 = pow(dist2,2) + pow(EPS,2);
                double pref2 = -G/pow(preff2,1.5)*galaxy.stars[j].mass;
                galaxy2.stars[i].acceleration.x += pref2 * dx2;
                galaxy2.stars[i].acceleration.y += pref2 * dy2;
                galaxy2.stars[i].acceleration.z += pref2 * dz2;
                continue;
            }
            double dx = galaxy2.stars[i].position.x - galaxy2.stars[j].position.x;
            double dy = galaxy2.stars[i].position.y - galaxy2.stars[j].position.y;
            double dz = galaxy2.stars[i].position.z - galaxy2.stars[j].position.z;
            double dist = sqrt(dx * dx + dy * dy + dz * dz);
            double preff = pow(dist,2) + pow(EPS,2);
            double pref = -G/pow(preff,1.5)*galaxy2.stars[j].mass;
            galaxy2.stars[i].acceleration.x += pref * dx;
            galaxy2.stars[i].acceleration.y += pref * dy;
            galaxy2.stars[i].acceleration.z += pref * dz;

            //s druhou galaxiou
            double dx2 = galaxy2.stars[i].position.x - galaxy.stars[j].position.x;
            double dy2 = galaxy2.stars[i].position.y - galaxy.stars[j].position.y;
            double dz2 = galaxy2.stars[i].position.z - galaxy.stars[j].position.z;
            double dist2 = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
            double preff2 = pow(dist2,2) + pow(EPS,2);
            double pref2 = -G/pow(preff2,1.5)*galaxy.stars[j].mass;
            galaxy2.stars[i].acceleration.x += pref2 * dx2;
            galaxy2.stars[i].acceleration.y += pref2 * dy2;
            galaxy2.stars[i].acceleration.z += pref2 * dz2;

        }
    }

}

void cal_velocity(int i){
    double G = 6.6742367e-11; // m^3.kg^-1.s^-2
    double orbital_vel;
    double distance;
    VECTOR v,norm_v;
//    double dx = galaxy.stars[i].position.x-galaxy.center.x;
//    double dy = galaxy.stars[i].position.y-galaxy.center.y;
//    double dz = galaxy.stars[i].position.z-galaxy.center.z;

    v.x = galaxy.center.x - galaxy.stars[i].position.x;
    v.y = galaxy.center.y - galaxy.stars[i].position.y;
    v.z = galaxy.center.z - galaxy.stars[i].position.z;

    double dis = v.x*v.x+v.y*v.y+v.z*v.z;
    distance = sqrt(dis);
    orbital_vel = G*galaxy.mass/distance/distance;
    norm_v = norm_vector(v);
    galaxy.stars[i].velocity.x += (norm_v.x*orbital_vel)*time_step;
    galaxy.stars[i].velocity.y += (norm_v.y*orbital_vel)*time_step;
    galaxy.stars[i].velocity.z += (norm_v.z*orbital_vel)*time_step;

}


