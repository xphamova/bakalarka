#include <GL/glut.h>
#include <GLFW/glfw3.h>
#include "galaxy.c"
#include "time.h"
#include "stdio.h"

void reshape(int, int);

void myInit();

void myDraw();

void gravity_calculate_acceleration();

void gravity_calculate_acceleration2();

void integrator_leapfrog();

void integrator_leapfrog2();


#define num_star 1000

double half_time_step;
double time_step;

GALAXY galaxy, galaxy2;

int main(int argc, char *argv[]) {

    srand((unsigned int) time(NULL));
    float height_frequency = rand_float_from_to(0, 1.5f);
    float height_magnitude = rand_float_from_to(0, 0.7f);

    //prva galaxia
    VECTOR galaxy_center, velocity;
    galaxy_center.x = 2;
    galaxy_center.y = 1;
    galaxy_center.z = 1;

    velocity.x = -15 * 1e5;
    velocity.y = 0;
    velocity.z = 0;

    galaxy = create_galaxy(height_magnitude, height_frequency, num_star, galaxy_center, velocity);

    //druha galaxia
    VECTOR galaxy_center_2, velocity2;
    galaxy_center_2.x = -2;
    galaxy_center_2.y = 1;
    galaxy_center_2.z = 1;

    velocity2.x = 15 * 1e5;
    velocity2.y = 0;
    velocity2.z = 0;

    galaxy2 = create_galaxy(height_magnitude, height_frequency, num_star, galaxy_center_2, velocity2);

    time_step = 1.0/32.0;
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

    gluLookAt(0.0f * 1e7, 10.0f * 1e7, -0.0f * 1e7, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
    glPointSize(1.0);
    glBegin(GL_POINTS);


     integrator_leapfrog();
//    integrator_leapfrog2();

    glColor3f(1.0, 1.0, 1.0);
    for (int i = 0; i < num_star; i++) {
        glVertex3f(galaxy.stars[i].position.x, galaxy.stars[i].position.y, galaxy.stars[i].position.z);
    }
    glColor3f(1.0, 0.2, 0.2);
    for (int i = 0; i < num_star; i++) {
        glVertex3f(galaxy2.stars[i].position.x, galaxy2.stars[i].position.y, galaxy2.stars[i].position.z);
    }

    glEnd();
    glutSwapBuffers();
}


void integrator_leapfrog() {
    for (int i = 0; i < num_star; i++) {
        //prva
        galaxy.stars[i].position.x += half_time_step * galaxy.stars[i].velocity.x;
        galaxy.stars[i].position.y += half_time_step * galaxy.stars[i].velocity.y;
        galaxy.stars[i].position.z += half_time_step * galaxy.stars[i].velocity.z;

        //druha
        galaxy2.stars[i].position.x += half_time_step * galaxy2.stars[i].velocity.x;
        galaxy2.stars[i].position.y += half_time_step * galaxy2.stars[i].velocity.y;
        galaxy2.stars[i].position.z += half_time_step * galaxy2.stars[i].velocity.z;
    }

    gravity_calculate_acceleration();
    gravity_calculate_acceleration2();
    for (int i = 0; i < num_star; i++) {
        //prva

        galaxy.stars[i].velocity.x += time_step * galaxy.stars[i].acceleration.x;
        galaxy.stars[i].velocity.y += time_step * galaxy.stars[i].acceleration.y;
        galaxy.stars[i].velocity.z += time_step * galaxy.stars[i].acceleration.z;

        galaxy.stars[i].position.x += half_time_step * galaxy.stars[i].velocity.x;
        galaxy.stars[i].position.y += half_time_step * galaxy.stars[i].velocity.y;
        galaxy.stars[i].position.z += half_time_step * galaxy.stars[i].velocity.z;

        //druha
        galaxy2.stars[i].velocity.x += time_step * galaxy2.stars[i].acceleration.x;
        galaxy2.stars[i].velocity.y += time_step * galaxy2.stars[i].acceleration.y;
        galaxy2.stars[i].velocity.z += time_step * galaxy2.stars[i].acceleration.z;

        galaxy2.stars[i].position.x += half_time_step * galaxy2.stars[i].velocity.x;
        galaxy2.stars[i].position.y += half_time_step * galaxy2.stars[i].velocity.y;
        galaxy2.stars[i].position.z += half_time_step * galaxy2.stars[i].velocity.z;

    }

}
//iny sposob ratania
void integrator_leapfrog2(){

    for(int i=0; i<num_star; i++){
        galaxy.stars[i].velocity.x += half_time_step * galaxy.stars[i].acceleration.x;
        galaxy.stars[i].velocity.y += half_time_step * galaxy.stars[i].acceleration.y;
        galaxy.stars[i].velocity.z += half_time_step * galaxy.stars[i].acceleration.z;

        galaxy2.stars[i].velocity.x += half_time_step * galaxy2.stars[i].acceleration.x;
        galaxy2.stars[i].velocity.y += half_time_step * galaxy2.stars[i].acceleration.y;
        galaxy2.stars[i].velocity.z += half_time_step * galaxy2.stars[i].acceleration.z;
    }
    for(int i=0; i<num_star; i++){
        galaxy.stars[i].position.x += time_step * galaxy.stars[i].velocity.x;
        galaxy.stars[i].position.y += time_step * galaxy.stars[i].velocity.y;
        galaxy.stars[i].position.z += time_step * galaxy.stars[i].velocity.z;

        galaxy2.stars[i].position.x += time_step * galaxy2.stars[i].velocity.x;
        galaxy2.stars[i].position.y += time_step * galaxy2.stars[i].velocity.y;
        galaxy2.stars[i].position.z += time_step * galaxy2.stars[i].velocity.z;
    }
    gravity_calculate_acceleration();
    gravity_calculate_acceleration2();

    for(int i=0; i<num_star; i++){
        galaxy.stars[i].velocity.x += half_time_step * galaxy.stars[i].acceleration.x;
        galaxy.stars[i].velocity.y += half_time_step * galaxy.stars[i].acceleration.y;
        galaxy.stars[i].velocity.z += half_time_step * galaxy.stars[i].acceleration.z;

        galaxy2.stars[i].velocity.x += half_time_step * galaxy2.stars[i].acceleration.x;
        galaxy2.stars[i].velocity.y += half_time_step * galaxy2.stars[i].acceleration.y;
        galaxy2.stars[i].velocity.z += half_time_step * galaxy2.stars[i].acceleration.z;
    }

}

//vypocet pre prvu
void gravity_calculate_acceleration() {
    double G = 6.6742367e-11; // m^3.kg^-1.s^-2

    for (int i = 0; i < num_star; i++) {

        double sum_x = 0, sum_y = 0, sum_z = 0;
        for (int j = 0; j < num_star; j++) {
            if (j == i) {

                double dx2 = galaxy.stars[i].position.x - galaxy2.stars[j].position.x;
                double dy2 = galaxy.stars[i].position.y - galaxy2.stars[j].position.y;
                double dz2 = galaxy.stars[i].position.z - galaxy2.stars[j].position.z;
                double r2 = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
                double a2 = -G*galaxy2.stars[j].mass/pow(r2,2);
                sum_x += a2/r2*dx2;
                sum_y += a2/r2*dy2;
                sum_z += a2/r2*dz2;
                continue;
            }
            //vypocet vzdielonosti medzi hviezdami vramci svojej galaxie
            double dx = galaxy.stars[i].position.x - galaxy.stars[j].position.x;
            double dy = galaxy.stars[i].position.y - galaxy.stars[j].position.y;
            double dz = galaxy.stars[i].position.z - galaxy.stars[j].position.z;
            double r = sqrt(dx * dx + dy * dy + dz * dz);
            double a = -G*galaxy.stars[j].mass/pow(r,2);
            sum_x += a/r*dx;
            sum_y += a/r*dy;
            sum_z += a/r*dz;

            //s druhou galaxiou
            double dx2 = galaxy.stars[i].position.x - galaxy2.stars[j].position.x;
            double dy2 = galaxy.stars[i].position.y - galaxy2.stars[j].position.y;
            double dz2 = galaxy.stars[i].position.z - galaxy2.stars[j].position.z;
            double r2 = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
            double a2 = -G*galaxy2.stars[j].mass/pow(r2,2);
            sum_x += a2/r2*dx2;
            sum_y += a2/r2*dy2;
            sum_z += a2/r2*dz2;

        }
        galaxy.stars[i].acceleration.x = sum_x;
        galaxy.stars[i].acceleration.y = sum_y;
        galaxy.stars[i].acceleration.z = sum_z;
    }

}

//vypocet pre druhu
void gravity_calculate_acceleration2() {
    double G = 6.6742367e-11; // m^3.kg^-1.s^-2

    for (int i = 0; i < num_star; i++) {
        double sum_x = 0, sum_y = 0, sum_z = 0;
        for (int j = 0; j < num_star; j++) {
            if (j == i) {
                double dx2 = galaxy2.stars[i].position.x - galaxy.stars[j].position.x;
                double dy2 = galaxy2.stars[i].position.y - galaxy.stars[j].position.y;
                double dz2 = galaxy2.stars[i].position.z - galaxy.stars[j].position.z;
                double r2 = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
                double a2 = -G*galaxy.stars[j].mass/pow(r2,2);
                sum_x += a2/r2*dx2;
                sum_y += a2/r2*dy2;
                sum_z += a2/r2*dz2;
                continue;
            }
            double dx = galaxy2.stars[i].position.x - galaxy2.stars[j].position.x;
            double dy = galaxy2.stars[i].position.y - galaxy2.stars[j].position.y;
            double dz = galaxy2.stars[i].position.z - galaxy2.stars[j].position.z;
            double r = sqrt(dx * dx + dy * dy + dz * dz);
            double a = -G*galaxy2.stars[j].mass/pow(r,2);
            sum_x += a/r*dx;
            sum_y += a/r*dy;
            sum_z += a/r*dz;

            //s druhou galaxiou
            double dx2 = galaxy2.stars[i].position.x - galaxy.stars[j].position.x;
            double dy2 = galaxy2.stars[i].position.y - galaxy.stars[j].position.y;
            double dz2 = galaxy2.stars[i].position.z - galaxy.stars[j].position.z;
            double r2 = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
            double a2 = -G*galaxy.stars[j].mass/pow(r2,2);
            sum_x += a2/r2*dx2;
            sum_y += a2/r2*dy2;
            sum_z += a2/r2*dz2;

        }
        galaxy2.stars[i].acceleration.x = sum_x;
        galaxy2.stars[i].acceleration.y = sum_y;
        galaxy2.stars[i].acceleration.z = sum_z;

    }

}