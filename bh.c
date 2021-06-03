#include "galaxy.c"
#include "stdio.h"
//ok
typedef struct {

    unsigned int bodies;

    VECTOR vector_low,vector_mid,vector_top;

    VECTOR position_star;

    struct OCTNODE *children[8];

    void *point_to;

    double node_mass;

}OCTNODE;
//ok
typedef struct {
    struct OCTNODE *root_node;
    int node;
}BARNESHUT;
//ok
typedef struct {
    double mass;
    VECTOR COM;
}BH_NODE;

int sub_insert(OCTNODE *node,VECTOR position,double ,void *);
//ok
struct OCTNODE * create_node(double x1, double y1, double z1, double x2, double y2, double z2){
    OCTNODE *node = malloc(sizeof(OCTNODE));
    if (!node)
        return NULL;
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

    node->position_star.x = 0;
    node->position_star.y = 0;
    node->position_star.z = 0;

    node->children[0] = NULL;
    node->children[1] = NULL;
    node->children[2] = NULL;
    node->children[3] = NULL;
    node->children[4] = NULL;
    node->children[5] = NULL;
    node->children[6] = NULL;
    node->children[7] = NULL;

    node->bodies = 0;
    node->point_to = NULL;
    return (struct OCTNODE *) node;
}
//ok
int insert(OCTNODE *node,VECTOR position_star,double mass,void *usr_val){
   if(!node)
       return 0;
    if (node->bodies == 0){
        node->position_star = position_star;
        node->node_mass = mass;
        node->point_to = usr_val;
    }
    else{
        if (node->bodies == 1) {
            sub_insert(node, node->position_star,mass,node->point_to);
            node->point_to = NULL;
        }else sub_insert(node,position_star,mass,usr_val);
    }
    node->bodies++;
     return node->bodies;
}
//ok
int sub_insert(OCTNODE *node,VECTOR position,double mass,void *usr_val){
    //nadelenie uzla na 8 pod uzlov, zistujeme v kotrom pod uzly sa nachadza bod
    if (!node)
        return 0;
    double min_x, min_y, min_z, max_x, max_y, max_z;
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
   // insert((OCTNODE *) node->children[sub], position, mass, point_to);
      return insert((OCTNODE *) node->children[sub], position,mass,usr_val);
}
//ok
BARNESHUT* BarnesHut_creat(double min_x,double min_y,double min_z,double max_x,double max_y,double max_z){
    BARNESHUT *barneshut = malloc(sizeof(BARNESHUT));
    if(!barneshut)
        return NULL;
    barneshut->root_node = create_node(min_x,min_y,min_z,max_x,max_y,max_z);
   if(!(barneshut->root_node)){
       free(barneshut);
       return NULL;
   }
    barneshut->node = 0;
    return barneshut;
}
//ok
int BARNESHUT_add(BARNESHUT *barneshut,VECTOR position,double mass){
    if (!barneshut)
        return 0;
    BH_NODE  *bh = malloc(sizeof (BH_NODE));
    if(!bh)
        return 0;
    bh->mass = mass;
    bh->COM.x = position.x;
    bh->COM.y = position.y;
    bh->COM.z = position.z;
    insert((OCTNODE *) barneshut->root_node, position,mass,bh);
    return 1;
}
//ok
void Barneshut_cal_tree(OCTNODE* node){
    if(!node)
        return;
    if(node->bodies == 1)
        return ;
    else{
        node->point_to = malloc(sizeof (BH_NODE));
        BH_NODE *Bh_node = (BH_NODE*)(node->point_to);
        Bh_node->mass = 0;Bh_node->COM.x = 0;Bh_node->COM.y = 0;Bh_node->COM.z = 0;
        for (int i=0; i<8 ; i++){
            if(!node->children[i])
                continue;
            Barneshut_cal_tree((OCTNODE *) node->children[i]);
            OCTNODE *child_node = (OCTNODE *) node->children[i];
            BH_NODE *child_bh_node = (BH_NODE*) child_node->point_to;
            double child_mass = child_bh_node->mass;
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

int calculate_force(OCTNODE *node,STAR *star){
    if (!node)
        return 0;
    double G = 6.6742367e-11; // m^3.kg^-1.s^-2
    double EPS = 3e3;
    star->force.x = 0;star->force.y = 0;star->force.z = 0;

    BH_NODE bhNode = * (BH_NODE*)(node->point_to);

    double div_x = (bhNode.COM.x-star->position.x);
    double div_y = (bhNode.COM.y-star->position.y);
    double div_z = (bhNode.COM.z-star->position.z);
    double radius_ = sqrt(pow(div_x,2)+pow(div_y,2)+pow(div_z,2));
    if (radius_ == 0)
        return 1;
    double width = ((node->vector_top.x-node->vector_low.x) + (node->vector_top.y-node->vector_low.y) + (node->vector_top.z-node->vector_low.z))/3;

    if (width/radius_ < 0.5){
       double denom = pow(radius_,2) + pow(EPS,2);
        star->force.x = (G * bhNode.mass * star->mass * (bhNode.COM.x-star->position.x))/pow(denom,1.5) * 5.0f;
        star->force.y = (G * bhNode.mass * star->mass * (bhNode.COM.y-star->position.y))/pow(denom,1.5) * 5.0f;
        star->force.z = (G * bhNode.mass * star->mass * (bhNode.COM.z-star->position.z))/pow(denom,1.5) * 5.0f;
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
    return 1;
}

//free funkcie

void free_barneshut_tree(OCTNODE *node){
    if(!node)
        return;
    free(node->point_to);
    free_barneshut_tree((OCTNODE *) node->children[0]);
    free_barneshut_tree((OCTNODE *) node->children[1]);
    free_barneshut_tree((OCTNODE *) node->children[2]);
    free_barneshut_tree((OCTNODE *) node->children[3]);
    free_barneshut_tree((OCTNODE *) node->children[4]);
    free_barneshut_tree((OCTNODE *) node->children[5]);
    free_barneshut_tree((OCTNODE *) node->children[6]);
    free_barneshut_tree((OCTNODE *) node->children[7]);

}
void free_node(OCTNODE *node) {
    if (!node) return;
    node->point_to = NULL;
    free_node((OCTNODE *) node->children[0]);
    free_node((OCTNODE *) node->children[1]);
    free_node((OCTNODE *) node->children[2]);
    free_node((OCTNODE *) node->children[3]);
    free_node((OCTNODE *) node->children[4]);
    free_node((OCTNODE *) node->children[5]);
    free_node((OCTNODE *) node->children[6]);
    free_node((OCTNODE *) node->children[7]);
    free(node);
}
void free_barneshut_tree_bh(BARNESHUT *bh){

    free_barneshut_tree_bh((BARNESHUT *) bh->root_node);
    free_node((OCTNODE *) bh->root_node);
    free(bh);
}



