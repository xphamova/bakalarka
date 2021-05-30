#include "galaxy.c"
#include "stdio.h"

typedef struct {

    unsigned int bodies;

    VECTOR vector_low,vector_mid,vector_top;

    VECTOR position_star;

    struct OCTNODE *children[8];

    void *usr_val;

    float node_mass;

}OCTNODE;

typedef struct {
    struct OCTNODE *root_node;
    int node;
}BARNESHUT;

typedef struct {
    float mass;
    VECTOR COM;
}BH_NODE;

void sub_insert(OCTNODE *node,VECTOR position,float ,void *);

struct OCTNODE * create_node(float x1, float y1, float z1, float x2, float y2, float z2){
    OCTNODE *node = malloc(sizeof(OCTNODE));
    // ohranicenie priestoru
    node->vector_top.x = (x1 < x2 ? x2 : x1);
    node->vector_low.x = (x1 < x2 ? x1 : x2);
    node->vector_mid.x = (x1 + x2 )/2;

    node->vector_top.y = (y1 < y2 ? y2 : y1);
    node->vector_low.y = (y1 < y2 ? y1 : y2);
    node->vector_mid.y = (y1 + y2 )/2;

    node->vector_top.z = (z1 < z2 ? z2 : z1);
    node->vector_low.z = (z1 < z2 ? z1 : z2);
    node->vector_mid.z = (z1 + z2 )/2;

    node->children[0] = NULL;
    node->children[1] = NULL;
    node->children[2] = NULL;
    node->children[3] = NULL;
    node->children[4] = NULL;
    node->children[5] = NULL;
    node->children[6] = NULL;
    node->children[7] = NULL;

    node->bodies = 0;

    return (struct OCTNODE *) node;
}

void insert(OCTNODE *node,VECTOR position_star,float mass,void *usr_val){
    if (node->bodies == 0){
        node->position_star = position_star;
        node->node_mass = mass;
        node->usr_val = usr_val;
    }
    else{
        if (node->bodies == 1) {
            sub_insert(node, node->position_star,mass,node->usr_val);
            node->usr_val = NULL;
        }else sub_insert(node,position_star,mass,usr_val);
    }
    node->bodies++;
    // return node->bodies;
}
void sub_insert(OCTNODE *node,VECTOR position,float mass,void *usr_val){
    //nadelenie uzla na 8 pod uzlov, zistujeme v kotrom pod uzly sa nachadza bod
    float min_x, min_y, min_z, max_x, max_y, max_z;
    int sub = 0;
    if(position.x > node->vector_mid.x){
        sub += 1;
        min_x = node->vector_mid.x;
        max_x = node->vector_top.x;
    } else {
        min_x = node->vector_low.x;
        max_x = node->vector_mid.x;
    }

    if (position.y > node->vector_mid.y){
        sub +=2;
        min_y = node->vector_mid.y;
        max_y = node->vector_top.y;
    } else {
        min_y = node->vector_low.y;
        max_y = node->vector_mid.y;
    }

    if (position.z > node->vector_mid.z){
        sub +=4;
        min_z = node->vector_mid.z;
        max_z = node->vector_top.z;
    } else {
        min_z = node->vector_low.z;
        max_z = node->vector_mid.z;
    }
    // ak neexistuje pod uzol, vytvori ho na zaklade toho kde sa hviezda nachadza
    if(!(node->children[sub]))
        (node->children[sub]) = create_node(min_x,min_y,min_z,max_x,max_y,max_z);

    //vratime uzol
    insert(node->children[sub],position,mass,usr_val);
    //  return insert((OCTNODE *) node->children[sub], position);
}

BARNESHUT* BarnesHut_creat(float min_x,float min_y,float min_z,float max_x,float max_y,float max_z){
    BARNESHUT *barneshut = malloc(sizeof (BARNESHUT));
    barneshut->root_node = create_node(min_x,min_y,min_z,max_x,max_y,max_z);
    barneshut->node = 0;
    return barneshut;
}

int BARNESHUT_add(BARNESHUT *barneshut,VECTOR position,float mass){
    BH_NODE  *bh = malloc(sizeof (BH_NODE));
    bh->mass = mass;
    bh->COM.x = position.x;
    bh->COM.y = position.y;
    bh->COM.z = position.z;
    insert((OCTNODE *) barneshut->root_node, position,mass,bh);
    return 1;
}

void Barneshut_cal_tree(OCTNODE* node){
    if(node->bodies == 1)
        return ;
    else{
        node->usr_val = malloc(sizeof (BH_NODE));
        BH_NODE *Bh_node = (BH_NODE*)(node->usr_val);
        Bh_node->mass = 0;
        Bh_node->COM.x = 0;
        Bh_node->COM.y = 0;
        Bh_node->COM.z = 0;
        for (int i=0; i<8 ; i++){
            //ak dieta neexistuje,pokracuj
            if(!node->children[i])
                continue;
            Barneshut_cal_tree((OCTNODE *) node->children[i]);
            OCTNODE *child_node = node->children[i];
            BH_NODE *child_bh_node = (BH_NODE*) child_node->usr_val;
            float child_mass = child_bh_node->mass;
            Bh_node->mass += child_mass;
            Bh_node->COM.x += child_mass*child_bh_node->COM.x;
            Bh_node->COM.y += child_mass*child_bh_node->COM.y;
            Bh_node->COM.z += child_mass*child_bh_node->COM.z;
        }
        Bh_node->COM.x /= Bh_node->mass;
        Bh_node->COM.y /= Bh_node->mass;
        Bh_node->COM.z /= Bh_node->mass;
    }

}

void calculate_force(OCTNODE *node,STAR *star){
    double G = 6.6742367e-11; // m^3.kg^-1.s^-2
    double EPS = 3e4;
    star->force.x = 0;
    star->force.y = 0;
    star->force.z = 0;

    //nacitame hodnoty uzla
    BH_NODE bhNode = * (BH_NODE*)(node->usr_val);
    //vypocet vzdialenosti medzi COM a bodom
    float div_x = (bhNode.COM.x-star->position.x);
    float div_y = (bhNode.COM.y-star->position.y);
    float div_z = (bhNode.COM.z-star->position.z);
    float radius_ = sqrtf(powf(div_x,2)+powf(div_y,2)+powf(div_z,2));
    //  float radius = sqrtf(powf(bhNode.COM.x,2)+powf(bhNode.COM.y,2)+powf(bhNode.COM.z,2));
    // ak sa radius = 0 ?
    // vypocet sirky uzla kvoli 3d vyratame priemer
    float width = ((node->vector_top.x-node->vector_low.x) + (node->vector_top.y-node->vector_low.y) + (node->vector_top.z-node->vector_low.z))/3;
    // prahova hodnota, vseobecne pouzivana 0.5
    if (width/radius_ < 0.5){

        //doplnenie
       // float radius_over_3 = powf(radius_,3);
       double denom = pow(radius_,2) + pow(EPS,2);
        star->force.x = (G * bhNode.mass * star->mass * (bhNode.COM.x-star->position.x))/pow(denom,1.5);
        star->force.y = (G * bhNode.mass * star->mass * (bhNode.COM.y-star->position.y))/pow(denom,1.5);
        star->force.z = (G * bhNode.mass * star->mass * (bhNode.COM.z-star->position.z))/pow(denom,1.5);
    } else {
        for (int i=0; i < 8; i++){
            STAR child_force;
            child_force.force.x = 0; child_force.force.y = 0; child_force.force.z = 0; child_force.mass = star->mass;
            child_force.position = star->position;
            if(node->children[i]){
                calculate_force((OCTNODE *) node->children[i], &child_force);
                star->force.x += child_force.force.x;
                star->force.y += child_force.force.y;
                star->force.z += child_force.force.z;
            }
        }
    }
}

//free funkcie

void free_barneshut_tree(){

}

void free_node(OCTNODE *node) {
    if (!node) return;
    node->usr_val = NULL;
    free_node((OCTNODE *) node->children[0]);
    free_node(node->children[1]);
    free_node(node->children[2]);
    free_node(node->children[3]);
    free_node(node->children[4]);
    free_node(node->children[5]);
    free_node(node->children[6]);
    free_node(node->children[7]);
    free(node);
}