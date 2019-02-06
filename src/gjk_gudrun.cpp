/*
 * Gudrun Thorkelsdottir
 * 27/06/2018
 * practice collision detection using GJK algorithm
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

//#define STANDALONE

//using namespace std;

using std::cout;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::string;

struct point {
  double x, y, z;
};

struct vector {
  double x, y, z;
};

int sLen = 12;
int tLen = 12;
point* start;               // = new point[sLen];
point* test;                // = new point[tLen];
point* simplex = new point[4]; //for a max of 4 points (tetrahedron)
int simS = 0;
vector D;       //direction vector
point o;        //origin


 void input(){
  string ignore, s, t;

  ifstream in("collision.xyz");                               //ONE FILE
  getline(in, ignore);
  getline(in, ignore);
  sLen = 12, tLen = 12;
  start = new point[sLen];
  test = new point[tLen];
  for(int i = 0; i < sLen; i++)
    in >> ignore >> start[i].x >> start[i].y >> start[i].z;
  for(int i = 0; i < tLen; i++)
    in >> ignore >> test[i].x >> test[i].y >> test[i].z;
  in.close();

  // ifstream inS("test1.xyz");                                  //TWO FILES
  // getline(inS, s);
  // istringstream holdS(s);
  // holdS >> sLen;
  // start = new point[sLen];
  // getline(inS, ignore);
  // for(int i = 0; i < sLen; i++)
  //   inS >> ignore >> start[i].x >> start[i].y >> start[i].z;
  // inS.close();

  // ifstream inT("test2.xyz");
  // getline(inT, t);
  // istringstream holdT(t);
  // holdT >> tLen;
  // test = new point[tLen];
  // getline(inT, ignore);
  // for(int i = 0; i < tLen; i++)
  //   inT >> ignore >> test[i].x >> test[i].y >> test[i].z;
  // inT.close();
 };

void initialize(){

#ifdef STANDALONE
  input();
#endif

  D.x = 1;       //setting direction vector to X
  D.y = 0;
  D.z = 0;

  o.x = 0;     //setting origin
  o.y = 0;
  o.z = 0;
};

double dot(point p, vector v){                        //dot product of point and vector
  return (p.x * v.x) + (p.y * v.y) + (p.z * v.z);
};

double dot(vector one, vector two){                   //dot product of two vectors
  return (one.x * two.x) + (one.y * two.y) + (one.z * two.z);
};

vector cross(vector one, vector two){                 //cross product of two vectors
  vector rtn;
  rtn.x = ((one.y*two.z) - (one.z*two.y));
  rtn.y = ((one.z*two.x) - (one.x*two.z));
  rtn.z = ((one.x*two.y) - (one.y*two.x));
  return rtn;
};
vector neg(vector v){                                //returns negative of a vector
  vector rtn;
  rtn.x = 0, rtn.y = 0, rtn.z = 0;
  if(v.x != 0)
    rtn.x = -(v.x);
  if(v.y != 0)
    rtn.y = -(v.y);
  if(v.z != 0)
    rtn.z = -(v.z);
  return rtn;
};

vector sub(point one, point two){          //returns the vector between two points
  vector rtn;
  rtn.x = two.x - one.x;
  rtn.y = two.y = one.y;
  rtn.z = two.z - one.z;
  return rtn;
};

void writeToFile(){                                   //writes the data to xyz files for visualization
  ofstream o;
  o.open("CollisionPoints.xyz");
  o << sLen  + tLen<< endl;
  o << "//comment line" << endl;
  for(int i = 0; i < sLen; i++)
    o << "oxygen " << start[i].x << " " << start[i].y << " " << start[i].z << endl;
  for(int i = 0; i < tLen; i++)
    o << "carbon " << test[i].x << " " << test[i].y << " " << test[i].z << endl;
  o.close();

};

point support(){                      //returns the points furthest along in direction D on the M.Diff
  point max;
  double maxS = dot(start[0], D);
  double maxT = dot(test[0], neg(D));
  int indexS = 0, indexT = 0;
  for(int i = 1; i < sLen; i++){
    if(dot(start[i], D) > maxS){
      maxS = dot(start[i], D);
      indexS = i;
    }
  }
  for(int i = 1; i < tLen; i++){
    if(dot(test[i], neg(D)) > maxT){
      maxT = dot(test[i], neg(D));
      indexT = i;
    }
  }
  max.x = start[indexS].x - test[indexT].x;
  max.y = start[indexS].y - test[indexT].y;
  max.z = start[indexS].z - test[indexT].z;
  return max;
};

bool checkSimplex(point p){               //if simplex is a tetrahedron that encloses the origin, returns true
                                   //else, modifies simplex, search direction, returns false
  switch(simS){
  case 2:{ //simplex is a line
    vector po = sub(p, o);
    D = neg(cross(cross(sub(p, simplex[0]), po), sub(p, simplex[0])));  //D = PX x PO x PX  //NEGATIVE?????
    return false;
  }
  case 3:{ //simplex is a triangle
    vector po = sub(p, o);
    vector pxy = cross(sub(p, simplex[0]), sub(p, simplex[1]));  //PXY = PX x PY
    if(dot(cross(pxy, sub(p, simplex[1])), po) > 0){  //if (PXY x PY) T PO > 0
      simplex[0] = simplex[1];
      simplex[1] = p;
      simplex[2] = o;
      D = cross(cross(sub(p, simplex[1]), po), sub(p, simplex[1])); //D = PY x PO x PY
      simS--;
    }//end big if

    else if(dot(cross(sub(p, simplex[0]), pxy), po) > 0){  //(PX x PXY) T PO > 0
      simplex[1] = p;
      simplex[2] = o;
      D = cross(cross(sub(p, simplex[0]), po), sub(p, simplex[0])); //D = PX x PO x PX
      simS--;
    }
    else{         //origin is inside triangle
      if(dot(pxy, po) > 0)      //origin is above triangle
	D = pxy;
      else                //origin is below triangle
    	D = neg(pxy);
    }
    return false;
  }
  case 4:{ // simplex is a tetrahedron
    vector po = sub(p, o);
    //vector pyx = cross(sub(p, simplex[1]), sub(p, simplex[0]));  //PXY = PX x PY
    vector pyx = cross(sub(p, simplex[0]), sub(p, simplex[1]));  //PXY = PX x PY
    vector pzy = cross(sub(p, simplex[2]), sub(p, simplex[1]));  //PXZ = PX x PZ
    vector pxz = cross(sub(p, simplex[0]), sub(p, simplex[2]));  //PYZ = PY x PZ

    double k, l;
    k = dot(p, pyx);
    l = dot(simplex[2] , pyx);
    if(k*l > 0){
      simplex[2] = p;
      simplex[3] = o;
      simS--;
    if((k < l && l < 0) || (k > l && l > 0))
    	D = neg(pyx);
    else
    	D = pyx;
      return false;
    }

    k = dot(p, pzy);
    l = dot(simplex[0], pzy);
    if(k*l > 0){
      simplex[0] = simplex[1];
      simplex[1] = simplex[2];
      simplex[2] = p;
      simplex[3] = o;
      simS--;
      if((k < l && l < 0) || (k > l && l > 0))
    	D = neg(pzy);
      else
    	D = pzy;
      return false;
    }

    k = dot(p, pxz);
    l = dot(simplex[1], pxz);
    if(k*l > 0){
      simplex[1] = simplex[2];
      simplex[2] = p;
      simplex[3] = o;
      simS--;
      if((k < l && l < 0) || (k > l && l > 0))
    	D = neg(pxz);
      else
    	D = pxz;
      return false;
    }
    return true;



    // if(dot(pyx, po) > 0){       //PYX T PO > 0 //trinagular case with pyx
    //   point x = simplex[0];
    //   simplex[0] = simplex[1];
    //   simplex[1] = x;
    //   simplex[2] = p;
    //   simplex[3] = o;
    //   D = pyx;
    //   simS --;
    //   return false;
    // }
    // else if(dot(pzy, po) > 0) {  //PZY T PO > 0 //triangular case with pzy

    //   simplex[0] = simplex[2];
    //   simplex[2] = p;
    //   simplex[3] = o;
    //   D = pzy;
    //   simS --;
    //   return false;
    // }
    // else if(dot(pxz, po) > 0) {  //PXZ T PO > 0 //triangular case with pxz
    //   simplex[1] = simplex[2];
    //   simplex[2] = p;
    //   simplex[3] = o;
    //   D = pxz;
    //   simS --;
    //   return false;
    // }
    // else  //origin inside tetrahedron
    //   return true;
  }
  default: 
    cout << "something went wrong with simplex" << endl;
  }//end switch

 return false;
}

bool algorithm(){           //returns true if there is an intersection, false if there is no intersection
  simplex[simS] = support();
  simS++;
  D = neg(D);
  int i = 0;
  while(1 == 1){
    point p = support();
    simplex[simS] = p;
    simS++;
    if(dot(p, D) < 0)      //CHECK THIS!!!!!!!!!!!!!!!
      return false;          //NO INTERSECTION
    else {
      if(checkSimplex(p))                //Simplex is checked -- if true then...
  	return true;                                                 //INTERSECTION CONFIRMED
    }
    i++;
  }
};

bool gjk_algorithm_gudrun(double* r1, int n1, double* r2, int n2){
  simS = 0;
  start = new point[sLen];
  test = new point[tLen];
  for(int s = 0; s < n1; s++){
    start[s].x = r1[3*s + 0];
    start[s].y = r1[3*s + 1];
    start[s].z = r1[3*s + 2];
  }
  for(int t = 0; t < n2; t++){
    test[t].x = r2[3*t + 0];
    test[t].y = r2[3*t + 1];
    test[t].z = r2[3*t + 2];
  }
  initialize();
  bool rtn = algorithm();
#ifdef STANDALONE
  writeToFile();
#endif
  delete start;
  delete test;
  return rtn;
};

#ifdef STANDALONE
int main(){
  initialize();
  if(algorithm())
    cout << "COLLISION DETECTED" << endl;
  else
    cout << "NO COLLISION DETECTED" << endl;
  writeToFile();

  delete start;
  delete test;
  return 0;
};
#endif
