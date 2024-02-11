#include "phylib.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>

phylib_object *phylib_new_still_ball(unsigned char number, phylib_coord *pos){

    phylib_object *newObj = (phylib_object*)malloc(sizeof(phylib_object));   

    if (newObj == NULL){   //success check
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }

    //constructors
    newObj -> type = PHYLIB_STILL_BALL;           
    newObj -> obj.still_ball.number = number;     
    newObj -> obj.still_ball.pos = *pos;

    //return pointer
    return newObj;

}

phylib_object *phylib_new_rolling_ball(unsigned char number, phylib_coord *pos, phylib_coord *vel, phylib_coord *acc){


    phylib_object *newObj = (phylib_object*)malloc(sizeof(phylib_object));   

    if (newObj == NULL){   //success check
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }

    //constructors
    newObj -> type = PHYLIB_ROLLING_BALL;
    newObj -> obj.rolling_ball.number = number;    
    newObj -> obj.rolling_ball.pos = *pos;
    newObj -> obj.rolling_ball.vel = *vel;
    newObj -> obj.rolling_ball.acc = *acc;

    return newObj;
}   

phylib_object *phylib_new_hole(phylib_coord *pos){

    phylib_object *newObj = (phylib_object*)malloc(sizeof(phylib_object));   

    if (newObj == NULL){   //success check
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }

    //constructors
    newObj -> type = PHYLIB_HOLE;
    newObj -> obj.hole.pos = *pos;


    return newObj;

}

phylib_object *phylib_new_hcushion(double y){

    phylib_object *newObj = (phylib_object*)malloc(sizeof(phylib_object));   

    if (newObj == NULL){   //success check
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }

    //constructors
    newObj -> type = PHYLIB_HCUSHION;
    newObj -> obj.hcushion.y = y;


    return newObj;
}

phylib_object *phylib_new_vcushion(double x){

    phylib_object *newObj = (phylib_object*)malloc(sizeof(phylib_object));   

    if (newObj == NULL){   //success check
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }

    //constructors
    newObj -> type = PHYLIB_VCUSHION;
    newObj -> obj.vcushion.x = x;

    return newObj;
}

phylib_table *phylib_new_table(void){

    phylib_table *newTable = (phylib_table*)malloc(sizeof(phylib_table));   

    if (newTable == NULL){   //success check
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }

    //constructors
    newTable->time = 0.0;

     // Assign values to array elements for cushions/sides of table
    newTable->object[0] = phylib_new_hcushion(0.0);  // Horizontal cushion at y=0.0
    newTable->object[1] = phylib_new_hcushion(PHYLIB_TABLE_LENGTH);  // Horizontal cushion at y=0.0
    newTable->object[2] = phylib_new_vcushion(0.0);  // Horizontal cushion at y=0.0
    newTable->object[3] = phylib_new_vcushion(PHYLIB_TABLE_WIDTH);  // Horizontal cushion at y=0.0

    //assign values to array elements for 6 holes around the table
    phylib_coord *newPos = (phylib_coord*)malloc(sizeof(phylib_coord));   //create new pos object
 
    //change coordinates, then malloc a hole and put it into the array 
    newPos -> x = 0.0;
    newPos -> y = 0.0;
    newTable -> object[4] = phylib_new_hole(newPos);

    newPos -> x = 0.0; 
    newPos -> y = PHYLIB_TABLE_LENGTH / 2.0;
    newTable -> object[5] = phylib_new_hole(newPos);

    newPos -> x = 0.0;
    newPos -> y =  PHYLIB_TABLE_LENGTH;
    newTable -> object[6] = phylib_new_hole(newPos);

    newPos -> x = PHYLIB_TABLE_WIDTH;
    newPos -> y = 0.0;
    newTable -> object[7] = phylib_new_hole(newPos);

    newPos -> x = PHYLIB_TABLE_WIDTH;
    newPos -> y = PHYLIB_TABLE_LENGTH / 2.0;
    newTable -> object[8] = phylib_new_hole(newPos);

    newPos -> x = PHYLIB_TABLE_WIDTH;
    newPos -> y = PHYLIB_TABLE_LENGTH;
    newTable -> object[9] = phylib_new_hole(newPos);

    //everything else null
    for(int i = 10; i < PHYLIB_MAX_OBJECTS; i++){
        newTable -> object[i] = NULL;
    }

    return newTable;

}

void phylib_copy_object(phylib_object **dest, phylib_object **src){


  if (*src == NULL) {
        *dest = NULL;
    } else {
        *dest = (phylib_object *)malloc(sizeof(phylib_object));  //make new object to store value
        if (*dest == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            return;
        }
        memcpy(*dest, *src, sizeof(phylib_object));  
    }

}

phylib_table *phylib_copy_table(phylib_table *table){

    phylib_table *newTable = (phylib_table*)malloc(sizeof(phylib_table));   

    if (newTable == NULL){   //success check
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }
    //transfer time
    newTable -> time = table -> time;

    //copy objects over
    for(int i = 0; i < PHYLIB_MAX_OBJECTS; i++){
        newTable -> object[i] = table -> object [i];
    }


    return newTable;

}


void phylib_add_object(phylib_table *table, phylib_object *object){

    for(int i = 0; i < PHYLIB_MAX_OBJECTS; i++){
        if(table -> object[i] == NULL){     //look for first open spot
            table-> object[i] = object;
            return;
        }
        }    
    }    


void phylib_free_table(phylib_table *table) {
    // Free each object in the table
    for (int i = 0; i < PHYLIB_MAX_OBJECTS; i++) {
        if (table->object[i] != NULL) {
            table->object[i] = NULL;
            free(table->object[i]);
        }
    }
    // Free the table itself
    free(table);
}



phylib_coord phylib_sub(phylib_coord c1, phylib_coord c2){

    phylib_coord newCoord = {0.0, 0.0};

    newCoord.x = c1.x - c2.x;
    newCoord.y = c1.y - c2.y;

    return newCoord;
}

double phylib_length(phylib_coord c){

    double res_squared = (c.x * c.x) + (c.y * c.y);
    double res = sqrt(res_squared);

    return res;
}

double phylib_dot_product(phylib_coord a, phylib_coord b){

    double resDot = (a.x * b.x) + (a.y * b.y);

    return resDot;
}

double phylib_distance(phylib_object *obj1, phylib_object *obj2){
    if(obj1 == NULL || obj1 -> type != PHYLIB_ROLLING_BALL){         //has to be rolling ball
        return 1.0;
    }

    double dx, dy, distance;


    if(obj2 != NULL){

        //calculate distances between each object, for each rolling ball
        if(obj2 -> type == PHYLIB_STILL_BALL){             
            dx = (obj2->obj.still_ball.pos.x) - (obj1->obj.rolling_ball.pos.x);        
            dy = (obj2->obj.still_ball.pos.y) - (obj1->obj.rolling_ball.pos.y);
            distance = sqrt(dx * dx + dy * dy) - PHYLIB_BALL_DIAMETER;
            return distance;
        }
        else if(obj2 -> type == PHYLIB_ROLLING_BALL){
            dx = (obj2->obj.rolling_ball.pos.x) - (obj1->obj.rolling_ball.pos.x);
            dy = (obj2->obj.rolling_ball.pos.y) - (obj1->obj.rolling_ball.pos.y);
            distance = sqrt(dx * dx + dy * dy) - PHYLIB_BALL_DIAMETER;
            return distance;
        }
        else if(obj2 -> type == PHYLIB_HOLE){
            dx = (obj2->obj.hole.pos.x) - (obj1->obj.rolling_ball.pos.x);
            dy = (obj2->obj.hole.pos.y) - (obj1->obj.rolling_ball.pos.y);
            distance = sqrt(dx * dx + dy * dy) - PHYLIB_HOLE_RADIUS;
            if(distance <= 0){
                (obj1 -> type = PHYLIB_STILL_BALL);
            }
            return distance;
        }
        else if(obj2 -> type == PHYLIB_HCUSHION){
            distance = fabs((obj2->obj.hcushion.y) - (obj1->obj.rolling_ball.pos.y));
            return distance;
        }
        else if(obj2 -> type == PHYLIB_VCUSHION){
            distance = fabs((obj2->obj.vcushion.x) - (obj1->obj.rolling_ball.pos.x));
            return distance;
        }
        else{        
            return 1.0;
        }
    }

    return 420.69;
}

void phylib_roll(phylib_object *new, phylib_object *old, double time){

    if(old -> type != PHYLIB_ROLLING_BALL || new -> type != PHYLIB_ROLLING_BALL){     //both have to be rolling balls
        return;
    }
 
        //update positions and velocities after set time
        new->obj.rolling_ball.pos.x = old->obj.rolling_ball.pos.x + (old->obj.rolling_ball.vel.x * time) + (0.5*(old->obj.rolling_ball.acc.x)* time * time);
        new->obj.rolling_ball.pos.y = old->obj.rolling_ball.pos.y + (old->obj.rolling_ball.vel.y * time) + (0.5*(old->obj.rolling_ball.acc.y)* time * time);
    
        new->obj.rolling_ball.vel.x = old->obj.rolling_ball.vel.x + (old->obj.rolling_ball.acc.x * time);
        new->obj.rolling_ball.vel.y = old->obj.rolling_ball.vel.y + (old->obj.rolling_ball.acc.y * time);

        //stop balls if wall hit
        if(old->obj.rolling_ball.vel.x * new->obj.rolling_ball.vel.x < 0){
            new->obj.rolling_ball.vel.x = 0;
            new->obj.rolling_ball.acc.x = 0;
        }
        if(old->obj.rolling_ball.vel.y * new->obj.rolling_ball.vel.y < 0){
            new->obj.rolling_ball.vel.y = 0;
            new->obj.rolling_ball.acc.y = 0;
        }     
    }


unsigned char phylib_stopped(phylib_object *object){

    //has to exist and be rolling ball
    if(object == NULL || object -> type != PHYLIB_ROLLING_BALL){
        return 0;
    }

    double speedSquared = ((object->obj.rolling_ball.vel.x * object->obj.rolling_ball.vel.x) + (object->obj.rolling_ball.vel.y * object->obj.rolling_ball.vel.y));
    
    //speed vector
    double speed = sqrt(speedSquared);

    //make ball still ball if below threshold
    if(speed < PHYLIB_VEL_EPSILON && object != NULL){
        object->type = PHYLIB_STILL_BALL;
        return 1;
    }


    return 0;
}


void phylib_bounce(phylib_object **a, phylib_object **b){
    phylib_object *obj_a = *a;
    phylib_object *obj_b = *b;

    //still becomes rolling
    if(obj_b->type == PHYLIB_STILL_BALL){
        obj_b->type = PHYLIB_ROLLING_BALL;
    }

    if(obj_b-> type == PHYLIB_HCUSHION){
        obj_a->obj.rolling_ball.acc.y *= -1;
        obj_a->obj.rolling_ball.vel.y *= -1;
    }
    else if(obj_b->type == PHYLIB_VCUSHION){
        obj_a->obj.rolling_ball.acc.x *= -1;
        obj_a->obj.rolling_ball.vel.x *= -1;
    }
    else if(obj_b->type == PHYLIB_HOLE){
        *a = NULL;
    }
    else if(obj_b->type == PHYLIB_ROLLING_BALL){
        
        obj_b->type = PHYLIB_ROLLING_BALL;
        //physics accoridng to given equations on assignment outline
        double r_ab_x = obj_a->obj.rolling_ball.pos.x - obj_b->obj.rolling_ball.pos.x;
        double r_ab_y = obj_a->obj.rolling_ball.pos.y - obj_b->obj.rolling_ball.pos.y;

        double v_rel_x = obj_a->obj.rolling_ball.vel.x - obj_b->obj.rolling_ball.vel.x;
        double v_rel_y = obj_a->obj.rolling_ball.vel.y - obj_b->obj.rolling_ball.vel.y;

        double length_r_ab = sqrt(r_ab_x * r_ab_x + r_ab_y * r_ab_y);

        double n_x = r_ab_x / length_r_ab;
        double n_y = r_ab_y / length_r_ab;

        double v_rel_n = v_rel_x * n_x + v_rel_y * n_y;

        obj_a->obj.rolling_ball.vel.x -= v_rel_n * n_x;
        obj_a->obj.rolling_ball.vel.y -= v_rel_n * n_y;

        obj_b->obj.rolling_ball.vel.x += v_rel_n * n_x;
        obj_b->obj.rolling_ball.vel.y += v_rel_n * n_y;

        double speed_a = sqrt(obj_a->obj.rolling_ball.vel.x * obj_a->obj.rolling_ball.vel.x + obj_a->obj.rolling_ball.vel.y * obj_a->obj.rolling_ball.vel.y);
        double speed_b = sqrt(obj_b->obj.rolling_ball.vel.x * obj_b->obj.rolling_ball.vel.x + obj_b->obj.rolling_ball.vel.y * obj_b->obj.rolling_ball.vel.y);

        if (speed_a > PHYLIB_VEL_EPSILON) {
          obj_a->obj.rolling_ball.acc.x = -obj_a->obj.rolling_ball.vel.x / speed_a * PHYLIB_DRAG;
          obj_a->obj.rolling_ball.acc.y = -obj_a->obj.rolling_ball.vel.y / speed_a * PHYLIB_DRAG;
        }

        if (speed_b > PHYLIB_VEL_EPSILON) {
          obj_b->obj.rolling_ball.acc.x = -obj_b->obj.rolling_ball.vel.x / speed_b * PHYLIB_DRAG;
          obj_b->obj.rolling_ball.acc.y = -obj_b->obj.rolling_ball.vel.y / speed_b * PHYLIB_DRAG;
        }
    }
}


unsigned char phylib_rolling(phylib_table *t){

    //# of rolling balls on table
    int rollingBalls = 0;

    for(int i = 0; i < PHYLIB_MAX_OBJECTS; i++){
        if((t->object[i] != NULL) && (t->object[i]->type == PHYLIB_ROLLING_BALL)){   
            rollingBalls++;
        }
    }
    return rollingBalls;
}

phylib_table *phylib_segment(phylib_table *table){
    if(phylib_rolling(table) == 0){
        
        return NULL;
    }
    //simulate balls as they roll
    phylib_table *copy_table = phylib_copy_table(table);

    double time = table->time;
    
    
    while(time <= PHYLIB_MAX_TIME){
        
        for(int i = 0; i < PHYLIB_MAX_OBJECTS; i++){
            
            phylib_object *new_obj = copy_table->object[i];
            
            //simulate rolling balls
            if(new_obj != NULL && new_obj->type == PHYLIB_ROLLING_BALL){
                phylib_roll(new_obj, table->object[i], PHYLIB_SIM_RATE);
            
           
                for(int j = 0; j < PHYLIB_MAX_OBJECTS; j++){
                        if(i != j){
                            //check for bounce between 2 balls
                            if(phylib_distance(copy_table->object[i], copy_table->object[j]) <= 0.0){
                                phylib_bounce(&new_obj, &copy_table->object[j]);

                                if(copy_table->object[j]->type == PHYLIB_HOLE){
                                    copy_table->object[i] = NULL;
                                }

                                return copy_table;
                            }
                        }
                }
            
            //return if ball stopped
            if(new_obj != NULL){
                if(phylib_stopped(new_obj)){
                    return copy_table;
                }
                
            }
            }
            
        }
        
        time += PHYLIB_SIM_RATE;
        copy_table->time = time;
    }
    
    return copy_table;
        
}

void phylib_print_object( phylib_object *object )
{
  if (object==NULL)
  {
    printf( "NULL;\n" );
    return;
  }

  switch (object->type)
  {
    case PHYLIB_STILL_BALL:
      printf( "STILL_BALL (%d,%6.1lf,%6.1lf)\n",
	      object->obj.still_ball.number,
	      object->obj.still_ball.pos.x,
	      object->obj.still_ball.pos.y );
      break;

    case PHYLIB_ROLLING_BALL:
      printf( "ROLLING_BALL (%d,%6.1lf,%6.1lf,%6.1lf,%6.1lf,%6.1lf,%6.1lf)\n",
              object->obj.rolling_ball.number,
              object->obj.rolling_ball.pos.x,
              object->obj.rolling_ball.pos.y,
              object->obj.rolling_ball.vel.x,
              object->obj.rolling_ball.vel.y,
              object->obj.rolling_ball.acc.x,
              object->obj.rolling_ball.acc.y );
      break;

    case PHYLIB_HOLE:
      printf( "HOLE (%6.1lf,%6.1lf)\n",
	      object->obj.hole.pos.x,
	      object->obj.hole.pos.y );
      break;

    case PHYLIB_HCUSHION:
      printf( "HCUSHION (%6.1lf)\n",
	      object->obj.hcushion.y );
      break;

    case PHYLIB_VCUSHION:
      printf( "VCUSHION (%6.1lf)\n",
	      object->obj.vcushion.x );
      break;
  }
}

void phylib_print_table( phylib_table *table )
{
  if (!table)
  {
    printf( "NULL\n" );
    return ;
  }

  printf( "time = %6.1lf;\n", table->time );
  for ( int i=0; i<PHYLIB_MAX_OBJECTS; i++ )
  {
    printf( "  [%02d] = ", i );
    phylib_print_object( table->object[i] );
  }

}


int main( int argc, char **argv )
{
  phylib_coord pos, vel, acc;
  phylib_table *table;
  phylib_object *sb;
  phylib_object *rb;

  table = phylib_new_table();

  // create a still ball 1/4 of the way "down" the middle of the table,
  // shift it up, and to the left just a little bit
  pos.x = PHYLIB_TABLE_WIDTH / 2.0 
          - sqrt( PHYLIB_BALL_DIAMETER*PHYLIB_BALL_DIAMETER / 2.0 );
  pos.y = PHYLIB_TABLE_WIDTH / 2.0
          - sqrt( PHYLIB_BALL_DIAMETER*PHYLIB_BALL_DIAMETER / 2.0 );
  sb = phylib_new_still_ball( 1, &pos );

  // create a rolling ball 3/4 of the way "down the table,
  // rolling up along the centre
  pos.x = PHYLIB_TABLE_WIDTH / 2.0;
  pos.y = PHYLIB_TABLE_LENGTH - PHYLIB_TABLE_WIDTH / 2.0;
  vel.x = 0.0;
  vel.y = -1000.0; // 1m/s (medium speed)
  acc.x = 0.0;
  acc.y = 180.0;
  rb = phylib_new_rolling_ball( 0, &pos, &vel, &acc );

  phylib_add_object( table, sb );
  phylib_add_object( table, rb );

  phylib_print_table( table );

  do
  {
    phylib_table *new = phylib_segment( table );
    phylib_free_table( table );
    table = new;

    phylib_print_table( table );
  } while( table );
}
