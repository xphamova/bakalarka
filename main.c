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

void reshape(int, int);

void myInit();

void myDraw();

void gravity_calculate_acceleration(int,int);

void gravity_calculate_acceleration2(int,int);

void *calculate_galaxy2(void*);

void* calculate_galaxy1(void*);

void *bh_start();

void start_cal();

STAR update_position(STAR);

STAR acceleration(STAR);

STAR second_update(STAR);

void start_thread();

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
    galaxy_center.x = 0;
    galaxy_center.y = 0;
    galaxy_center.z = 0;

   //velocity.x = 15 * 1e5;
    velocity.x = 0;
    velocity.y = 0;
    velocity.z = 0;


    galaxy = create_galaxy(height_magnitude, height_frequency, num_star, galaxy_center, velocity);

    //druha galaxia
    VECTOR galaxy_center_2, velocity2;
    galaxy_center_2.x = 1.5;
    galaxy_center_2.y = 0;
    galaxy_center_2.z = 0;

//    velocity2.x = -15 * 1e5;
    velocity2.x = 0;
    velocity2.y = 0;
    velocity2.z = 0;

    galaxy2 = create_galaxy(height_magnitude, height_frequency, num_star, galaxy_center_2, velocity2);

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

   // start_cal();
   start_thread();
    glColor3f(1.0, 1.0, 1.0);
    for (int i = 0; i < num_star; i++) {
        glVertex3f(galaxy.stars[i].position.x, galaxy.stars[i].position.y, galaxy.stars[i].position.z);
    }
    glColor3f(1.0, 0.2, 0.2);
    for (int i = 0; i < num_star; i++) {
        glVertex3f(galaxy2.stars[i].position.x, galaxy2.stars[i].position.y, galaxy2.stars[i].position.z);
    }
    glVertex3f(0, 0, 0);
    glColor3f(1.0, 1.0, 1.0);
    glPointSize(5.0);
    glVertex3f(galaxy2.center.x, galaxy2.center.y, galaxy2.center.z);

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
    pthread_create(&tid3, NULL, calculate_galaxy2, (void *) send_range1);
    pthread_create(&tid4, NULL, calculate_galaxy2, (void *) send_range2);

    pthread_join(tid1, (void **) &receive_range1);
    pthread_join(tid2, (void **) &receive_range2);
    pthread_join(tid3, (void **) &receive_range3);
    pthread_join(tid4, (void **) &receive_range4);


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
    gravity_calculate_acceleration(start,end);

    for(int i = start; i<end; i++){

        galaxy.stars[i].velocity.x += time_step * galaxy.stars[i].acceleration.x;
        galaxy.stars[i].velocity.y += time_step * galaxy.stars[i].acceleration.y;
        galaxy.stars[i].velocity.z += time_step * galaxy.stars[i].acceleration.z;


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


    gravity_calculate_acceleration2(start,end);
    for (int i = start; i < end; i++) {
        galaxy2.stars[i].velocity.x += time_step * galaxy2.stars[i].acceleration.x;
        galaxy2.stars[i].velocity.y += time_step * galaxy2.stars[i].acceleration.y;
        galaxy2.stars[i].velocity.z += time_step * galaxy2.stars[i].acceleration.z;

        galaxy2.stars[i].position.x += half_time_step * galaxy2.stars[i].velocity.x;
        galaxy2.stars[i].position.y += half_time_step * galaxy2.stars[i].velocity.y;
        galaxy2.stars[i].position.z += half_time_step * galaxy2.stars[i].velocity.z;
    }
}

//vypocet pre prvu
void gravity_calculate_acceleration(int start, int end) {
    double G = 6.6742367e-11; // m^3.kg^-1.s^-2
    double EPS =3e4;

    for (int i = start; i < end; i++) {
        galaxy.stars[i].acceleration.x = 0;
        galaxy.stars[i].acceleration.y = 0;
        galaxy.stars[i].acceleration.z = 0;
        for (int j = 0; j < num_star; j++) {
            if (j == i) {
                double dx2 = galaxy.stars[i].position.x - galaxy2.stars[j].position.x;
                double dy2 = galaxy.stars[i].position.y - galaxy2.stars[j].position.y;
                double dz2 = galaxy.stars[i].position.z - galaxy2.stars[j].position.z;
                double dist2 = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
                double preff2 = pow(dist2,2) + pow(EPS,2);
                double pref2 = -G/pow(preff2,1.5)*galaxy2.stars[j].mass;
                galaxy.stars[i].acceleration.x += pref2 * dx2;
                galaxy.stars[i].acceleration.y += pref2 * dy2;
                galaxy.stars[i].acceleration.z += pref2 * dz2;
                continue;
            }
//            vypocet vzdielonosti medzi hviezdami vramci svojej galaxie
            double dx = galaxy.stars[i].position.x - galaxy.stars[j].position.x;
            double dy = galaxy.stars[i].position.y - galaxy.stars[j].position.y;
            double dz = galaxy.stars[i].position.z - galaxy.stars[j].position.z;

            double dist = sqrt(dx * dx + dy * dy+ dz * dz);
            double preff = pow(dist,2) + pow(EPS,2);
            double pref = -G/pow(preff,1.5)*galaxy.stars[j].mass;
            galaxy.stars[i].acceleration.x += pref * dx;
            galaxy.stars[i].acceleration.y += pref * dy;
            galaxy.stars[i].acceleration.z += pref * dz;

            //s druhou galaxiou
            double dx2 = galaxy.stars[i].position.x - galaxy2.stars[j].position.x;
            double dy2 = galaxy.stars[i].position.y - galaxy2.stars[j].position.y;
            double dz2 = galaxy.stars[i].position.z - galaxy2.stars[j].position.z;
            double dist2 = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
            double preff2 = pow(dist2,2) + pow(EPS,2);
            double pref2 = -G/pow(preff2,1.5)*galaxy2.stars[j].mass;
            galaxy.stars[i].acceleration.x += pref2 * dx2;
            galaxy.stars[i].acceleration.y += pref2 * dy2;
            galaxy.stars[i].acceleration.z += pref2 * dz2;

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
void start_thread(){

    //vytvorime strom pre vsetky hviezdy
    //   BARNESHUT *BH = BarnesHut_creat(-960,-540,-500,960,540,500);
    BARNESHUT *BH = BarnesHut_creat(-960*1e7,-540*1e7,-500*1e7,960*1e7,540*1e7,500*1e7);
    for (int i = 0;i<num_star;i++){
        BARNESHUT_add(BH,galaxy.stars[i].position,galaxy.stars[i].mass);
        BARNESHUT_add(BH,galaxy2.stars[i].position,galaxy2.stars[i].mass);
    }

    // vyratame COM kazdeho uzla
    Barneshut_cal_tree( BH->root_node);

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

    //podelienie hviezd vlaknam
    //vyriesit problem ak cislo nie je delitelne 4
    int number_of_star_in_thread = num_star/2;
    STAR range1[number_of_star_in_thread];
    STAR range2[number_of_star_in_thread],range3[number_of_star_in_thread],range4[number_of_star_in_thread];

    for (int i=0; i<number_of_star_in_thread;i++) {
        range1[i].mass = 0;
        range1[i].force.x = 0;
        range1[i].force.y = 0;
        range1[i].force.z = 0;
        range1[i].acceleration.x = 0;
        range1[i].acceleration.y = 0;
        range1[i].acceleration.z = 0;
        range1[i].velocity.x = 0;
        range1[i].velocity.y = 0;
        range1[i].velocity.z = 0;
        range1[i].position.x = 0;
        range1[i].position.y = 0;
        range1[i].position.z = 0;
    }
    for (int i=0; i<number_of_star_in_thread;i++) {
        range2[i].mass = 0;
        range2[i].force.x = 0;
        range2[i].force.y = 0;
        range2[i].force.z = 0;
        range2[i].acceleration.x = 0;
        range2[i].acceleration.y = 0;
        range2[i].acceleration.z = 0;
        range2[i].velocity.x = 0;
        range2[i].velocity.y = 0;
        range2[i].velocity.z = 0;
        range2[i].position.x = 0;
        range2[i].position.y = 0;
        range2[i].position.z = 0;
    }

    for (int i=0; i<number_of_star_in_thread;i++) {
        range4[i].mass = 0;
        range4[i].force.x = 0;
        range4[i].force.y = 0;
        range4[i].force.z = 0;
        range4[i].acceleration.x = 0;
        range4[i].acceleration.y = 0;
        range4[i].acceleration.z = 0;
        range4[i].velocity.x = 0;
        range4[i].velocity.y = 0;
        range4[i].velocity.z = 0;
        range4[i].position.x = 0;
        range4[i].position.y = 0;
        range4[i].position.z = 0;
    }
    int r1 = 0, r2=0, r3 = 0, r4=0;

    for(int i = 0; i<num_star; i++){
        if (i<number_of_star_in_thread){
            range1[r1] = galaxy.stars[i];
            r1++;
        }
        if (i<number_of_star_in_thread){
            range3[r3] = galaxy2.stars[i];
            r3++;
        }
        if (i>=number_of_star_in_thread ){
            range2[r2] = galaxy.stars[i];
            r2++;
        }
        if(i>=number_of_star_in_thread){
            range4[r4] = galaxy2.stars[i];
            r4++;
        }
    }

    void *send_range1,*send_range2,*send_range3,*send_range4;
    STAR *receive_range1,*receive_range2,*receive_range3,*receive_range4;
    send_range1 = &range1;
    send_range2 = &range2;
    send_range3 = &range3;
    send_range4 = &range4;
    pthread_t tid1, tid2,tid3,tid4;
    pthread_create(&tid1, NULL, bh_start, send_range1);//start thread
    pthread_create(&tid2, NULL, bh_start, (void *) send_range2);
    pthread_create(&tid3, NULL, bh_start, (void *) send_range3);
    pthread_create(&tid4, NULL, bh_start, (void *) send_range4);
    pthread_join(tid1, (void **) &receive_range1);
    pthread_join(tid2, (void **) &receive_range2);
    pthread_join(tid3, (void **) &receive_range3);
    pthread_join(tid4, (void **) &receive_range4);

    r1 = 0; r2=0; r3=0; r4=0;

    for(int i = 0; i<num_star; i++){
        if (i<number_of_star_in_thread){
            galaxy.stars[i]  = receive_range1[r1];
            r1++;
        }
        if (i<number_of_star_in_thread){
            galaxy2.stars[i] = receive_range3[r3];
            r3++;
        }
        if (i>= number_of_star_in_thread ){
            galaxy.stars[i] = receive_range2[r2];
            galaxy2.stars[i] = receive_range4[r2];
            r2++;
        }
        if(i>= number_of_star_in_thread){

        }
    }

    free_node((OCTNODE *) BH->root_node);
    free(BH);
    BH = NULL;
}

void *bh_start(void *input){

    STAR range[num_star/4];
    STAR *ran = (STAR *)input;
    int size = 0;

    for (int i = 0; i< num_star/4; i++){
        range[i] = *(ran + i);
        if(range[i].mass != 0){
            size++;
        }
    }

    //ziskanie stromu z pamate
    //shared memory
    key_t key = 26;
    BARNESHUT *BH;
    int shmid = shmget(key, sizeof(BARNESHUT), 0666);
    BARNESHUT *shm;
    shm = (BARNESHUT *) shmat(shmid, NULL, 0);
    BH = shm;



    for(int i = 0; i < num_star/4; i++){
        //kick : pos = pos + half_time_step * vel
        range[i] = update_position(range[i]);
        calculate_force( BH->root_node, &range[i]);
        //drift : acc[i] = range.force[i]/range.mass[i]
        range[i] = acceleration(range[i]);
        // kick : vel =vel + acc[i]* timestep; pos = half_time_step * vel
        range[i] = second_update(range[i]);
    }

    STAR* output;

    //create dynamic memory
    output = (STAR *) malloc(size * sizeof(STAR));
    if (output == NULL) {
        printf("Memory not allocated.\n");
        exit(0);
    }

    //write to dynamic memory
    for (int i = 0; i < size; i++) {
        output[i] = range[i];
    }

    /* detach from the segment: */
    if (shmdt(shm) == -1) {
        perror("shmdt");
        exit(1);
    }

    return (STAR*) output;

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

