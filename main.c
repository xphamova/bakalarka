#include <GL/glut.h>
#include <GLFW/glfw3.h>
//#include "galaxy.c"
#include "time.h"
#include "stdio.h"
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include "bh.c"
#include <sys/shm.h>
#include <signal.h>

void reshape(int, int);

void myInit();

void myDraw();

void bh_start();

STAR update_position(STAR);

STAR acceleration(STAR);

STAR second_update(STAR);

void start_thread();

void update_center_first();

void update_center_sec();

STAR accel_from_center(STAR);

void exitfunc(int);
#define num_star 5000

double half_time_step;
double time_step;

GALAXY galaxy, galaxy2;
int screen = 0;
int main(int argc, char *argv[]) {

    signal(SIGALRM, exitfunc);
    alarm(120);

    srand((unsigned int) time(NULL));
    float height_frequency = rand_float_from_to(0, 1.5f);
    float height_magnitude = rand_float_from_to(0, 0.7f);

    //prva galaxia
    VECTOR galaxy_center, velocity;
    galaxy_center.x = -0.5;
    galaxy_center.y = 0;
    galaxy_center.z = 0;

   velocity.x = 10 * 1e6;
  //  velocity.x = 0;
    velocity.y = 0;
    velocity.z = 0;


    galaxy = create_galaxy(height_magnitude, height_frequency, num_star, galaxy_center, velocity);

    //druha galaxia
    VECTOR galaxy_center_2, velocity2;
    galaxy_center_2.x = 0.5;
    galaxy_center_2.y = 0;
    galaxy_center_2.z = 0;

    velocity2.x = -10 * 1e6;
//    velocity2.x = 0;
    velocity2.y = 0;
    velocity2.z = 0;

    galaxy2 = create_galaxy(height_magnitude, height_frequency, num_star, galaxy_center_2, velocity2);

    time_step = 1.0/2048.0;
    half_time_step = 0.5 * time_step;

    /* Initialize window system */
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
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


   start_thread();
   update_center_first();
   update_center_sec();
    glColor3f(1.0, 1.0, 1.0);
    for (int i = 0; i < num_star; i++) {
        glVertex3f(galaxy.stars[i].position.x, galaxy.stars[i].position.y, galaxy.stars[i].position.z);
    }
    glColor3f(0.9, 0.9, 0.9);
    for (int i = 0; i < num_star; i++) {
        glVertex3f(galaxy2.stars[i].position.x, galaxy2.stars[i].position.y, galaxy2.stars[i].position.z);
    }
    glColor3f(0.9, 0.2, 0.2);
    glPointSize(15.0);

    glVertex3f(galaxy.center.x,galaxy.center.y, galaxy.center.z);
    glColor3f(0.9, 0.2, 0.2);
    glPointSize(15.0*1e7);
    glVertex3f(galaxy2.center.x, galaxy2.center.y, galaxy2.center.z);

    glEnd();
    glutSwapBuffers();

    screen+= 1;
    printf("%d\n",screen);
}

void start_thread(){

    //vytvorime strom pre vsetky hviezdy
    //   BARNESHUT *BH = BarnesHut_creat(-960,-540,-500,960,540,500);
    BARNESHUT *BH = BarnesHut_creat(-96*1e9,-54*1e9,-50*1e9,96*1e9,54*1e9,50*1e9);
    for (int i = 0;i<num_star;i++){
        BARNESHUT_add(BH,galaxy.stars[i].position,galaxy.stars[i].mass);
        BARNESHUT_add(BH,galaxy2.stars[i].position,galaxy2.stars[i].mass);
    }

    // vyratame COM kazdeho uzla
    Barneshut_cal_tree((OCTNODE *) BH->root_node);

    //zapiseme tieto veci do zdielanej pamate
    //shared memory
    key_t key = 26;
    //osetrenie na -1
    int shmid = shmget(key, sizeof(BARNESHUT), 0666 | IPC_CREAT);//identifikator
    //osetrenie na -1
    BARNESHUT *shm = (BARNESHUT *) shmat(shmid, NULL, 0);//prepajanie


    BARNESHUT *s;
    s = shm;
    *s = *BH;

    /* detach from the segment: */
    if (shmdt(shm) == -1) {
        perror("shmdt");
        exit(1);
    }

    bh_start();

    free_node((OCTNODE *) BH->root_node);
    free(BH);

    BH = NULL;
}

void bh_start(){

    int j = 0;
    STAR range[num_star*2];
    for (int i=0; i< num_star;i++){
        range[i]=galaxy.stars[i];
    }
    for (int i=num_star; i<num_star*2;i++){
        range[i] = galaxy2.stars[j];
        j++;
    }
    //ziskanie stromu z pamate
    //shared memory
    key_t key = 26;
    BARNESHUT *BH;
    int shmid = shmget(key, sizeof(BARNESHUT), 0666);
    BARNESHUT *shm;
    shm = (BARNESHUT *) shmat(shmid, NULL, 0);
    BH = shm;


    for(int i = 0; i < num_star*2; i++){
        //kick : pos = pos + half_time_step * vel
        range[i] = update_position(range[i]);
        calculate_force((OCTNODE *) BH->root_node, &range[i]);
        //drift : acc[i] = range.force[i]/range.mass[i]
        range[i] = acceleration(range[i]);
        range[i] = accel_from_center(range[i]);
        // kick : vel =vel + acc[i]* timestep; pos = half_time_step * vel
        range[i] = second_update(range[i]);
    }


    /* detach from the segment: */
    if (shmdt(shm) == -1) {
        perror("shmdt");
        exit(1);
    }
    int k=0;
    for(int i=0;i<num_star;i++){
        galaxy.stars[i]=range[i];
    }
    for(int i=num_star;i<num_star*2;i++){
        galaxy2.stars[k]=range[i];
        k++;
    }


}

//LEAPFROG------------------------------------------------------------------------------------------------


STAR update_position( STAR star){

    star.position.x += half_time_step*star.velocity.x;
    star.position.y += half_time_step*star.velocity.x;
    star.position.z += half_time_step*star.velocity.z;

    return star;
}

STAR acceleration(STAR star){
    star.acceleration.x = star.force.x/star.mass;
    star.acceleration.y = star.force.y/star.mass;
    star.acceleration.z = star.force.z/star.mass;
    return star;
}

STAR second_update(STAR star){

    star.velocity.x += time_step * star.acceleration.x;
    star.velocity.y += time_step * star.acceleration.y;
    star.velocity.z += time_step * star.acceleration.z;

    star.position.x += half_time_step * star.velocity.x;
    star.position.y += half_time_step * star.velocity.y;
    star.position.z += half_time_step * star.velocity.z;
    return star;

}

void update_center_first(){
    double G = 6.6742367e-11; // m^3.kg^-1.s^-2
    double EPS =3e4;
    galaxy.acceleration.x = 0;
    galaxy.acceleration.y = 0;
    galaxy.acceleration.z = 0;

    galaxy.center.x += half_time_step * galaxy.velocity.x;
    galaxy.center.y += half_time_step * galaxy.velocity.y;
    galaxy.center.z += half_time_step * galaxy.velocity.z;


    double dx = galaxy.center.x - galaxy2.center.x;
    double dy = galaxy.center.y - galaxy2.center.y;
    double dz = galaxy.center.z - galaxy2.center.z;
    double dist = sqrt(dx * dx + dy * dy + dz * dz);
    double preff = pow(dist,2) + pow(EPS,2);
    double pref = -G/pow(preff,1.5)*galaxy2.mass;
    galaxy.acceleration.x += pref * dx;
    galaxy.acceleration.y += pref * dy;
    galaxy.acceleration.z += pref * dz;


    galaxy.velocity.x += time_step * galaxy.acceleration.x;
    galaxy.velocity.y += time_step * galaxy.acceleration.y;
    galaxy.velocity.z += time_step * galaxy.acceleration.z;


    galaxy.center.x += half_time_step * galaxy.velocity.x;
    galaxy.center.y += half_time_step * galaxy.velocity.y;
    galaxy.center.z += half_time_step * galaxy.velocity.z;
}
void update_center_sec(){
    double G = 6.6742367e-11; // m^3.kg^-1.s^-2
    double EPS =3e4;
    galaxy2.acceleration.x = 0;
    galaxy2.acceleration.y = 0;
    galaxy2.acceleration.z = 0;

    galaxy2.center.x += half_time_step * galaxy2.velocity.x;
    galaxy2.center.y += half_time_step * galaxy2.velocity.y;
    galaxy2.center.z += half_time_step * galaxy2.velocity.z;


    double dx = galaxy2.center.x - galaxy.center.x;
    double dy = galaxy2.center.y - galaxy.center.y;
    double dz = galaxy2.center.z - galaxy.center.z;
    double dist = sqrt(dx * dx + dy * dy + dz * dz);
    double preff = pow(dist,2) + pow(EPS,2);
    double pref = -G/pow(preff,1.5)*galaxy.mass;
    galaxy2.acceleration.x += pref * dx;
    galaxy2.acceleration.y += pref * dy;
    galaxy2.acceleration.z += pref * dz;


    galaxy2.velocity.x += time_step * galaxy2.acceleration.x;
    galaxy2.velocity.y += time_step * galaxy2.acceleration.y;
    galaxy2.velocity.z += time_step * galaxy2.acceleration.z;


    galaxy2.center.x += half_time_step * galaxy2.velocity.x;
    galaxy2.center.y += half_time_step * galaxy2.velocity.y;
    galaxy2.center.z += half_time_step * galaxy2.velocity.z;
}

STAR accel_from_center(STAR star){
    double G = 6.6742367e-11; // m^3.kg^-1.s^-2
    double EPS =3e5;
    float num = 5.0f;
    double dx = star.position.x - galaxy.center.x;
    double dy = star.position.y - galaxy.center.y;
    double dz = star.position.z - galaxy.center.z;
    double dist = sqrt(dx * dx + dy * dy + dz * dz);
    double preff = pow(dist,2) + pow(EPS,2);
    double pref = -G/pow(preff,1.5)*galaxy.mass;
    star.acceleration.x += pref * dx * num;
    star.acceleration.y += pref * dy * num;
    star.acceleration.z += pref * dz * num;

    double dx2 = star.position.x - galaxy2.center.x;
    double dy2 = star.position.y - galaxy2.center.y;
    double dz2 = star.position.z - galaxy2.center.z;
    double dist2 = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
    double preff2 = pow(dist2,2) + pow(EPS,2);
    double pref2 = -G/pow(preff2,1.5)*galaxy2.mass;
    star.acceleration.x += pref2 * dx2 * num;
    star.acceleration.y += pref2 * dy2 * num;
    star.acceleration.z += pref2 * dz2 * num;
    return star;
}

void exitfunc(int sig)
{
    _exit(sig);
}