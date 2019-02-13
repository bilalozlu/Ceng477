#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <climits>

typedef unsigned char RGB[3];
using namespace parser;

parser::Scene scene;
Camera camera;
Vec3f e;
float distance;
float l;
float r;
float b;
float t;
int n_x;
int n_y;


// Sample usage for reading an XML scene file
enum ObjectType {
    S,
    T,
    M
};

// vector<unsigned char>* data;

float take_determinant (float a00,
                        float a01,
                        float a02,
                        float a10,
                        float a11,
                        float a12,
                        float a20,
                        float a21,
                        float a22) {

                        return a00*(a11*a22-a12*a21)+a10*(a02*a21-a01*a22)+a20*(a01*a12-a11*a02);

}

//(x,y,z)/s
Vec3f vector_division(Vec3f v,float s) {
  return {v.x/s, v.y/s, v.z/s };
}

//(x,y,z).(a,b,c) = (x*a,y*b,z*c)
Vec3f multiply_two_vectors(Vec3f v1,Vec3f v2) {
  return {v1.x*v2.x, v1.y*v2.y, v1.z*v2.z };
}

//(x,y,z)*s
Vec3f vector_mul(Vec3f v,float s) {
  return {v.x*s, v.y*s, v.z*s };
}


Vec3f vector_sum(Vec3f v1,Vec3f v2) {
  return {v1.x+v2.x, v1.y+v2.y, v1.z+v2.z };
}

Vec3f vector_subtraction(Vec3f v1,Vec3f v2) {
  return {v1.x-v2.x, v1.y-v2.y, v1.z-v2.z };
}

Vec3f cross_product(Vec3f v1, Vec3f v2) {
    return {v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x};
}

float dot_product(Vec3f v1, Vec3f v2) {
  return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}

float norm(Vec3f v) {
  return sqrt(dot_product(v,v));
}

// (v1.v2)/(|v1||v2|)d
float cosVal(Vec3f v1, Vec3f v2) {
  return fmax(dot_product(v1,v2)/(sqrt(dot_product(v1,v1)) * sqrt(dot_product(v2,v2))),0);
}

bool intersect_object_light(Vec3f p, Vec3f d, float distance) {

  for(Sphere sphere:scene.spheres) {
    Vec3f c = scene.vertex_data[sphere.center_vertex_id-1];
    float temp1 = (d.x * (p.x- c.x)) + (d.y * (p.y- c.y)) + (d.z * (p.z- c.z));
    float temp2 = dot_product(d,d);
    float temp3 = (p.x - c.x)*(p.x - c.x) + (p.y - c.y)*(p.y - c.y) + (p.z - c.z)*(p.z - c.z) - pow(sphere.radius,2);
    float determinant = pow(temp1,2) - temp2*temp3;
    if (determinant >= 0) {

      float ray_distance_to_light1 = (-temp1 - sqrt(determinant))/temp2;
      float ray_distance_to_light2 = (-temp1 + sqrt(determinant))/temp2;

      float temp4 = ray_distance_to_light1 < ray_distance_to_light2 ? ray_distance_to_light1:ray_distance_to_light2;

      if ((ray_distance_to_light1 < distance || ray_distance_to_light2 < distance ) &&  ray_distance_to_light1 > 0.01 && ray_distance_to_light2 > 0.01){
        return false;
      }
    }
  }

  for(Triangle triangle:scene.triangles) {

    Vec3f a = scene.vertex_data[triangle.indices.v0_id-1];
    Vec3f b = scene.vertex_data[triangle.indices.v1_id-1];
    Vec3f c = scene.vertex_data[triangle.indices.v2_id-1];

    float det_A = take_determinant(a.x-b.x,a.x-c.x,d.x,a.y-b.y,a.y-c.y,d.y,a.z-b.z,a.z-c.z,d.z);
    float beta = take_determinant(a.x-p.x,a.x-c.x,d.x,a.y-p.y,a.y-c.y,d.y,a.z-p.z,a.z-c.z,d.z)/det_A;
    float gama = take_determinant(a.x-b.x,a.x-p.x,d.x,a.y-b.y,a.y-p.y,d.y,a.z-b.z,a.z-p.z,d.z)/det_A;
    if(beta + gama <= 1 && 0 <= beta && 0<=gama) {
      float ray_distance_to_light = take_determinant(a.x-b.x, a.x-c.x, a.x-p.x, a.y-b.y, a.y-c.y, a.y-p.y, a.z-b.z, a.z-c.z, a.z-p.z)/det_A;
      if (ray_distance_to_light < distance && ray_distance_to_light>=0.01) {
        return false;
      }
    }
  }

  for(Mesh mesh:scene.meshes) {

    for(Face face:mesh.faces) {
      Vec3f a = scene.vertex_data[face.v0_id-1];
      Vec3f b = scene.vertex_data[face.v1_id-1];
      Vec3f c = scene.vertex_data[face.v2_id-1];

      float det_A = take_determinant(a.x-b.x,a.x-c.x,d.x,a.y-b.y,a.y-c.y,d.y,a.z-b.z,a.z-c.z,d.z);
      float beta = take_determinant(a.x-p.x,a.x-c.x,d.x,a.y-p.y,a.y-c.y,d.y,a.z-p.z,a.z-c.z,d.z)/det_A;
      float gama = take_determinant(a.x-b.x,a.x-p.x,d.x,a.y-b.y,a.y-p.y,d.y,a.z-b.z,a.z-p.z,d.z)/det_A;

      if(beta + gama <= 1 && 0 <= beta && 0<=gama) {
        float ray_distance_to_light = take_determinant(a.x-b.x, a.x-c.x, a.x-p.x, a.y-b.y, a.y-c.y, a.y-p.y, a.z-b.z, a.z-c.z, a.z-p.z)/det_A;

        if (ray_distance_to_light < distance && ray_distance_to_light>=0.01){
          return false;
        }
      }
    }
  }

  return true;

}

void find_triangle_normals(parser::Scene *scene) {
    int i=0;
    for(Triangle triangle:scene->triangles) {
      Vec3f a = scene->vertex_data[triangle.indices.v0_id-1];
      Vec3f b = scene->vertex_data[triangle.indices.v1_id-1];
      Vec3f c = scene->vertex_data[triangle.indices.v2_id-1];
      Vec3f n = cross_product(vector_subtraction(c,b), vector_subtraction(a,b));
      // std::cout<
      (scene->triangles)[i].normal_vector = vector_division(n,norm(n));

      i++;
    }
}

void find_mesh_normals(parser::Scene *scene) {
  int i=0;
  for(Mesh mesh:scene->meshes) {
    for(Face face:mesh.faces) {
      Vec3f a = scene->vertex_data[face.v0_id-1];
      Vec3f b = scene->vertex_data[face.v1_id-1];
      Vec3f c = scene->vertex_data[face.v2_id-1];
      Vec3f n = cross_product(vector_subtraction(c,b), vector_subtraction(a,b));

      (scene->meshes)[i].normal_vectors.push_back(vector_division(n,norm(n)));
    }
    i++;
  }

}
Vec3f find_sphere_normal(Vec3f intersected, Vec3f center) {
    Vec3f n = vector_subtraction(intersected, center);
    return vector_division(n,norm(n) );
}

//find intersecting material and minimum distance to it
Vec3i intersect(int &intersecting_material_id, Vec3f &intersecting_point,Vec3f &normal, Vec3f d, Vec3f e ,bool val) {

  float min = (float) INT_MAX;
  enum ObjectType objectType;
  Vec3f sphere_center;
  Vec3f triangle_normal;

  for(Sphere sphere:scene.spheres) {

    Vec3f c = scene.vertex_data[sphere.center_vertex_id -1 ];
    // find determinant: sqrt((d(o-c))^2 -d^2((o-c)^2 -R^2))
    //temp1 = (d(o-c))
    float temp1 = (d.x * (e.x- c.x)) + (d.y * (e.y- c.y)) + (d.z * (e.z- c.z));
    //temp2 = d.d
    float temp2 = d.x*d.x + d.y*d.y + d.z*d.z;
    //temp3 = (o-c)^2 -R^2
    float temp3 = (e.x - c.x)*(e.x - c.x) + (e.y - c.y)*(e.y - c.y) + (e.z - c.z)*(e.z - c.z) - pow(sphere.radius,2);
    float determinant = pow(temp1,2) - temp2*temp3;


    if(determinant >= 0) {
      float ray_distance_to_object1 = (-temp1 - sqrt(determinant))/temp2;
      float ray_distance_to_object2 = (-temp1 + sqrt(determinant))/temp2;

      float temp4 = ray_distance_to_object1 < ray_distance_to_object2 ? ray_distance_to_object1:ray_distance_to_object2;


      if(min > temp4 && temp4 >= 0.01) {
        min = temp4;
        intersecting_material_id = sphere.material_id;
        objectType = S;
        sphere_center = c;
      }
    }
  }

  for(Triangle triangle:scene.triangles) {

    Vec3f a = scene.vertex_data[triangle.indices.v0_id-1];
    Vec3f b = scene.vertex_data[triangle.indices.v1_id-1];
    Vec3f c = scene.vertex_data[triangle.indices.v2_id-1];

    float det_A = take_determinant(a.x-b.x,a.x-c.x,d.x,a.y-b.y,a.y-c.y,d.y,a.z-b.z,a.z-c.z,d.z);
    float beta = take_determinant(a.x-e.x,a.x-c.x,d.x,a.y-e.y,a.y-c.y,d.y,a.z-e.z,a.z-c.z,d.z)/det_A;
    float gama = take_determinant(a.x-b.x,a.x-e.x,d.x,a.y-b.y,a.y-e.y,d.y,a.z-b.z,a.z-e.z,d.z)/det_A;
    if(beta + gama <= 1.0 && 0.0 <= beta && 0.0<=gama) {

      float ray_distance_to_object = take_determinant(a.x-b.x, a.x-c.x, a.x-e.x, a.y-b.y, a.y-c.y, a.y-e.y, a.z-b.z, a.z-c.z, a.z-e.z)/det_A;


      if(min>ray_distance_to_object  && ray_distance_to_object >= 0.01) {
        min = ray_distance_to_object;
        intersecting_material_id = triangle.material_id;
        objectType = T;
        triangle_normal = triangle.normal_vector;

      }
    }
  }

  for(Mesh mesh:scene.meshes) {
    int i = 0;
    for(Face face:mesh.faces) {
      Vec3f a = scene.vertex_data[face.v0_id-1];
      Vec3f b = scene.vertex_data[face.v1_id-1];
      Vec3f c = scene.vertex_data[face.v2_id-1];

      float det_A = take_determinant(a.x-b.x,a.x-c.x,d.x,a.y-b.y,a.y-c.y,d.y,a.z-b.z,a.z-c.z,d.z);


        float beta = take_determinant(a.x-e.x,a.x-c.x,d.x,a.y-e.y,a.y-c.y,d.y,a.z-e.z,a.z-c.z,d.z)/det_A;
        float gama = take_determinant(a.x-b.x,a.x-e.x,d.x,a.y-b.y,a.y-e.y,d.y,a.z-b.z,a.z-e.z,d.z)/det_A;

        if(beta + gama <= 1.0 && 0.0 <= beta && 0.0<=gama) {
          float ray_distance_to_object = take_determinant(a.x-b.x, a.x-c.x, a.x-e.x, a.y-b.y, a.y-c.y, a.y-e.y, a.z-b.z, a.z-c.z, a.z-e.z)/det_A;

          //std::cout<<"ray_distance_to_object:"<<ray_distance_to_object<<std::endl;

          if(min > ray_distance_to_object && ray_distance_to_object >= 0.01) {
            min = ray_distance_to_object;
            intersecting_material_id = mesh.material_id;
            objectType = M;
            triangle_normal = (mesh.normal_vectors)[i];

          }
        }
        i++;



    }
  }

  if(intersecting_material_id == -1 ) {
    //set background color
    return {0,0,0};
  }
  else {

    intersecting_point.x = e.x + d.x*min;
    intersecting_point.y = e.y + d.y*min;
    intersecting_point.z = e.z + d.z*min;

    Material material = scene.materials[intersecting_material_id - 1 ];

    // total color variable
    Vec3i pixel_color = {0,0,0};

    // ambient shading
    Vec3f ambient_color_contribution = multiply_two_vectors(material.ambient,scene.ambient_light );
    Vec3f diffuse_color_contribution = {0,0,0};
    Vec3f specular_color_contribution = {0,0,0};

    //  point light sources seen by this point
    for(PointLight point_light:scene.point_lights) {

      Vec3f object_to_light_vector = vector_subtraction(point_light.position,intersecting_point);

      // object_to_light_vector = vector_mul(object_to_light_vector,scene.shadow_ray_epsilon);
      float light_object_distance = norm(object_to_light_vector);
      object_to_light_vector = vector_division(object_to_light_vector,norm(object_to_light_vector));

      intersecting_point = vector_subtraction(point_light.position,vector_mul(object_to_light_vector,light_object_distance + scene.shadow_ray_epsilon));


      if(intersect_object_light(intersecting_point, object_to_light_vector, light_object_distance) ) {

          // find half vector
          Vec3f temp = vector_sum(vector_mul(d,-1), object_to_light_vector);
          Vec3f half_vector = vector_division(temp,norm(temp));

          // specular_color_contribution = (I/d^2) * cos(theta)^R * (material_color)
          if(objectType == S) {
            normal = find_sphere_normal(intersecting_point, sphere_center);
          }
          else if(objectType == T || objectType == M) {
            normal = triangle_normal;
          }
          else {
            std::cout<<"error";
          }


          diffuse_color_contribution = vector_sum(diffuse_color_contribution, multiply_two_vectors (vector_division( point_light.intensity, pow(light_object_distance,2) ),
                                        vector_mul(material.diffuse,cosVal(normal,object_to_light_vector)) ));
          specular_color_contribution = vector_sum(specular_color_contribution, multiply_two_vectors(vector_division( point_light.intensity, pow(light_object_distance,2) ),
                                         vector_mul(material.specular, pow(cosVal(normal,half_vector),material.phong_exponent))));

        }
      }

      // cast to total color variable
      pixel_color.x = (int) (diffuse_color_contribution.x + ambient_color_contribution.x + specular_color_contribution.x + 0.5);
      pixel_color.y = (int) (diffuse_color_contribution.y + ambient_color_contribution.y + specular_color_contribution.y + 0.5);
      pixel_color.z = (int) (diffuse_color_contribution.z + ambient_color_contribution.z + specular_color_contribution.z + 0.5);

      return pixel_color;
    }
}

int main(int argc, char* argv[]){
  scene.loadFromXml(argv[1]);

  int s = scene.cameras.size();

  for(int i=0;i<s;i++) {
      // std::string image_name = camera.image_name;
      Camera camera = (scene.cameras)[i];
      find_triangle_normals(&scene);
      find_mesh_normals(&scene);

      Vec3f e = camera.position;
      distance = camera.near_distance;
      l = camera.near_plane.x;
      r = camera.near_plane.y;
      b = camera.near_plane.z;
      t = camera.near_plane.w;
      n_x = camera.image_width;
      n_y = camera.image_height;

      Vec3f v = camera.up;

      Vec3f w = vector_mul(camera.gaze,-1);
      Vec3f u ={v.y*w.z - v.z*w.y, v.z*w.x - v.x*w.z, v.x*w.y - v.y*w.x};

      unsigned char array[n_x*n_y*3];

      Vec3f m =  vector_subtraction(e,vector_mul(w,distance));
      Vec3f q = vector_sum(vector_sum(m,vector_mul(u,l)),vector_mul(v,t));

      for(int j=0;j<n_y;j++) {
        for(int i=0;i<n_x;i++) {

          float su = ((r - l) * (i+0.5)) /n_x;
          float sv = ((t - b) * (j+0.5)) /n_y;
          Vec3f s = vector_sum(q,vector_subtraction(vector_mul(u,su),vector_mul(v,sv)));

          Vec3f d = vector_subtraction(s,e);
          d = vector_division(d,norm(d));

          // std::cout<<"d: "<<d.x<<" "<<d.y<<" "<<d.z<<std::endl;

          int intersecting_material_id = -1;
          Vec3f intersecting_point = {-1,-1,-1};
          Vec3f normal;
          Vec3i pixel_color = intersect(intersecting_material_id,intersecting_point,normal, d,e,false);
          Vec3i mirror_color = {0,0,0};


          for(int depth = scene.max_recursion_depth;depth>0;depth--) {

            if(intersecting_material_id != -1) {
              // if(intersecting_material_id != 1)
              // std::cout<<"id: "<<intersecting_material_id<<std::endl;
              Material material = scene.materials[intersecting_material_id - 1 ];
              if(material.mirror.x > 0 || material.mirror.y >0 || material.mirror.z >0) {

                 // -wo + 2n(n.wo) , where w0 = -d
                d = vector_sum(d, vector_mul(normal,2*cosVal(normal, vector_mul(d,-1))));
                d = vector_division(d,norm(d));

                intersecting_material_id = -1;
                Vec3i intersect_color = intersect(intersecting_material_id,intersecting_point, normal, d, intersecting_point,true);
                mirror_color.x += intersect_color.x * material.mirror.x;
                mirror_color.y += intersect_color.y * material.mirror.y;
                mirror_color.z += intersect_color.z * material.mirror.z;

              }
              else{
                break;
              }
            }
            else {
              break;
            }

          }

          pixel_color.x += mirror_color.x;
          pixel_color.y += mirror_color.y;
          pixel_color.z += mirror_color.z;

          array[j*3*n_x + i*3] = pixel_color.x > 255?255:pixel_color.x;
          array[j*3*n_x + i*3+1] = pixel_color.y > 255?255:pixel_color.y;
          array[j*3*n_x + i*3+2] = pixel_color.z > 255?255:pixel_color.z;
        }
      }
        write_ppm(camera.image_name.c_str(),array,n_x,n_y);
  }


      return 0;
    }
