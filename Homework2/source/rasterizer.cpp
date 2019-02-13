#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "hw2_types.h"
#include "hw2_math_ops.h"
#include "hw2_file_ops.h"
#include <iostream>

#define PI 3.14159265


Camera cameras[100];
int numberOfCameras = 0;

Model models[1000];
int numberOfModels = 0;

Color colors[100000];
int numberOfColors = 0;

Translation translations[1000];
int numberOfTranslations = 0;

Rotation rotations[1000];
int numberOfRotations = 0;

Scaling scalings[1000];
int numberOfScalings = 0;

Vec3 vertices[100000];
int numberOfVertices = 0;

Color backgroundColor;

// backface culling setting, default disabled    // Vertex processing

int backfaceCullingSetting = 0;

Color **image;



/*
	Initializes image with background color
*/
void initializeImage(Camera cam) {
    int i, j;

    for (i = 0; i < cam.sizeX; i++)
        for (j = 0; j < cam.sizeY; j++) {
            image[i][j].r = backgroundColor.r;
            image[i][j].g = backgroundColor.g;
            image[i][j].b = backgroundColor.b;

        }
}

/*
	Transformations, culling, rasterization are done here.
	You can define helper functions inside this file (rasterizer.cpp) only.
	Using types in "hw2_types.h" and functions in "hw2_math_ops.cpp" will speed you up while working.
*/
double min3(double n1, double n2, double n3) {
        if (n2 >= n3 && n1 >= n3 ) {
                return n3;
        }
        else if (n1 >= n2 && n3 >= n2) {
                return n2;
        }
        else {
                return n1;
        }
}

double max3(double n1, double n2, double n3) {
        if (n1 >= n2 && n1 >= n3 ) {
                return n1;
        }
        else if (n2 >= n1 && n2 >= n3) {
                return n2;
        }
        else {
                return n3;
        }
}

double fxy(double x, double y, double x0, double x1, double y0, double y1){
        return x*(y0-y1) + y*(x1-x0) + x0*y1 -y0*x1;
}

void translation(const Translation &T, const Vec3 &in, Vec3 &out) {
        double r[4];
        double m1[4][4] = {{1, 0, 0, T.tx},
                          {0, 1, 0, T.ty},
                          {0, 0, 1, T.tz},
                          {0, 0, 0, 1}};
        double m2[4] = {in.x, in.y, in.z, 1};

        multiplyMatrixWithVec4d(r, m1, m2);

        out = {r[0], r[1], r[2]};
}

void scaling(const Scaling &S , const Vec3 &in , Vec3 &out) {
        double r[4];
        double m1[4][4] = {{S.sx, 0, 0, 0},
                          {0, S.sy, 0, 0},
                          {0, 0, S.sz, 0},
                          {0, 0, 0, 1}};
        double m2[4] = {in.x, in.y, in.z, 1};

        multiplyMatrixWithVec4d(r, m1, m2);

        out = {r[0], r[1], r[2]};
}

void rotation_transform(const Rotation &R, double out[4][4]) {
        double r1[4][4];
        double smallest = min3(R.ux, R.uy, R.uz);
        Vec3 V ;

        if(smallest == R.ux) {
                V.x = 0;
                V.y = R.uz;
                V.z = -R.uy;
        }
        else if(smallest == R.uy) {
                V.y = 0;
                V.x = R.uz;
                V.z = -R.ux;
        }
        else {
                V.z = 0;
                V.y = R.ux;
                V.x = -R.uy;
        }

        Vec3 W = crossProductVec3({R.ux, R.uy, R.uz}, V);

        double m1[4][4] = {{R.ux, R.uy, R.uz, 0},
                          {V.x, V.y, V.z, 0},
                          {W.x, W.y, W.z, 0},
                          {0, 0, 0, 1}};
        double m1_inverse[4][4] = {{R.ux, V.x, W.x, 0},
                                  {R.uy, V.y, W.y, 0},
                                  {R.uz, V.z, W.z, 0},
                                  {0, 0, 0, 1}};
        double rota1[4][4] = {{1, 0, 0, 0},
                             {0, cos(R.angle * PI / 180.0), -sin(R.angle * PI / 180.0), 0},
                             {0, sin(R.angle * PI / 180.0), cos(R.angle * PI / 180.0), 0},
                             {0, 0, 0, 1}};

        multiplyMatrixWithMatrix(r1, m1_inverse, rota1);
        multiplyMatrixWithMatrix(out, r1, m1);
}

void midpointLine(double x1, double x0, double y1, double y0, Vec3 color0, Vec3 color1, Camera cam ) {
        double alpha;
        double slope = (y1 - y0) / (x1 - x0);
        if (slope >= 1) {
                int x = x0;
                double d = 2*(x0 - x1) + y1 - y0;
                for (int y = y0; y <= y1; y++){
                        alpha = (x - x0) / (x1 - x0);

                        double cx = color0.x * (1 - alpha) + alpha * color1.x;
                        double cy = color0.y * (1 - alpha) + alpha * color1.y;
                        double cz = color0.z * (1 - alpha) + alpha * color1.z;

                        if (x > 0 && y > 0 && cam.sizeX > x && cam.sizeY > y){
                                image[x][y].r = (int)cx;
                                image[x][y].g = (int)cy;
                                image[x][y].b = (int)cz;
                        }

                        if (d < 0){
                                x++;
                                d += 2*(x0 - x1 + y1 - y0);
                        }
                        else{
                                d += 2*(x0 - x1);
                        }
                }
        }
        else if (slope > 0 ) {
                int y = y0;
                double d = 2*(y0 - y1) + x1 - x0;
                for (int x = x0; x <= x1; x++){
                        alpha = (x - x0) / (x1 - x0);

                        double cx = color0.x * (1 - alpha) + alpha * color1.x;
                        double cy = color0.y * (1 - alpha) + alpha * color1.y;
                        double cz = color0.z * (1 - alpha) + alpha * color1.z;

                        if (x > 0 && y > 0 && cam.sizeX > x && cam.sizeY > y){
                                image[x][y].r = (int)cx;
                                image[x][y].g = (int)cy;
                                image[x][y].b = (int)cz;
                        }

                        if (d < 0){
                                y++;
                                d += 2*(y0 - y1 + x1 - x0);
                        }
                        else{
                                d += 2*(y0 - y1);
                        }
                }
        }
        else if (slope >= -1){
                int y = y0;
                double d = 2*(y0 - y1) + x1 - x0;
                for (int x = x0; x <= x1; x++){
                        alpha = (x - x0) / (x1 - x0);

                        double cx = color0.x * (1 - alpha) + alpha * color1.x;
                        double cy = color0.y * (1 - alpha) + alpha * color1.y;
                        double cz = color0.z * (1 - alpha) + alpha * color1.z;

                        if (x > 0 && y > 0 && cam.sizeX > x && cam.sizeY > y){
                                image[x][y].r = (int)cx;
                                image[x][y].g = (int)cy;
                                image[x][y].b = (int)cz;
                        }

                        if (d < 0){
                                y--;
                                d += 2*(y1 - y0 + x1 - x0);
                        }
                        else{
                                d += 2*(y1 - y0);
                        }
                }
        }
        else {
                int x = x0;
                double d = 2*(x0 - x1) + y1 - y0;
                for (int y = y0; y >= y1; y--){
                        alpha = (x - x0) / (x1 - x0);

                        double cx = color0.x * (1 - alpha) + alpha * color1.x;
                        double cy = color0.y * (1 - alpha) + alpha * color1.y;
                        double cz = color0.z * (1 - alpha) + alpha * color1.z;

                        if (x > 0 && y > 0 && cam.sizeX > x && cam.sizeY > y){
                                image[x][y].r = (int)cx;
                                image[x][y].g = (int)cy;
                                image[x][y].b = (int)cz;
                        }

                        if (d < 0){
                                x++;
                                d += 2*(x0 - x1 + y0 - y1);
                        }
                        else{
                                d += 2*(x0 - x1);
                        }
                }
        }
}

void forwardRenderingPipeline(Camera cam) {
        bool isTransformed[numberOfVertices+1];
        Vec3 changedVertices[numberOfVertices+1];
        Vec3 changedVertices2[numberOfVertices+1];

        for (int a = 1; a <= numberOfVertices; a++) {
                changedVertices[a] = vertices[a];
        }

        //Modeling transformations
        for (int j = 0; j < numberOfModels; j++) {
                Model model = models[j];
                for (int i = 0; i < model.numberOfTransformations; i++) {
                        for (int k = 1; k <= numberOfVertices; k++) {
                                isTransformed[k] = false;
                        }

                        char trans_type = model.transformationTypes[i];
                        int trans_id = model.transformationIDs[i];

                        if (trans_type == 't') {
                                for (int k = 0; k < model.numberOfTriangles; k++) {
                                        Triangle cur_tri = model.triangles[k];
                                        for (int l = 0; l < 3; l++) {
                                                int cur_vertex_id = cur_tri.vertexIds[l];

                                                if (!isTransformed[cur_vertex_id]) {
                                                        Vec3 vertex = changedVertices[cur_vertex_id];
                                                        isTransformed[cur_vertex_id] = true;
                                                        translation(translations[trans_id], vertex, changedVertices[cur_vertex_id]);

                                                }
                                        }
                                }
                        }
                        else if (trans_type == 's') {
                                for (int k = 0; k < model.numberOfTriangles; k++) {
                                        Triangle cur_tri = model.triangles[k];
                                        for (int l = 0; l < 3; l++) {
                                                int cur_vertex_id = cur_tri.vertexIds[l];
                                                if (!isTransformed[cur_vertex_id]) {
                                                        Vec3 vertex = changedVertices[cur_vertex_id];
                                                        isTransformed[cur_vertex_id] = true;
                                                        scaling(scalings[trans_id], vertex, changedVertices[cur_vertex_id]);
                                                }
                                        }
                                }
                        }
                        else if (trans_type == 'r') {
                                for (int k = 0; k < model.numberOfTriangles; k++) {
                                        Triangle cur_tri = model.triangles[k];
                                        double theMatrix[4][4];
                                        rotation_transform(rotations[trans_id], theMatrix);
                                        for (int l = 0; l < 3; l++) {
                                                int cur_vertex_id = cur_tri.vertexIds[l];
                                                if (!isTransformed[cur_vertex_id]) {
                                                        double res[4];
                                                        Vec3 vertex = changedVertices[cur_vertex_id];
                                                        isTransformed[cur_vertex_id] = true;
                                                        double theVector[4] = {vertex.x, vertex.y, vertex.z, 1};
                                                        multiplyMatrixWithVec4d(res, theMatrix, theVector);
                                                        changedVertices[cur_vertex_id] = {res[0], res[1], res[2]};
                                                }
                                        }
                                }
                        }
                }
        }

        double Mcam[4][4];

        double R[4][4] = {{cam.u.x, cam.u.y, cam.u.z, 0},
                         {cam.v.x, cam.v.y, cam.v.z, 0},
                         {cam.w.x, cam.w.y, cam.w.z, 0},
                         {0, 0, 0, 1}};

        double E[4][4] = {{1, 0, 0, -cam.pos.x},
                         {0, 1, 0, -cam.pos.y},
                         {0, 0, 1, -cam.pos.z},
                         {0, 0, 0, 1}};

        double n = cam.n, f = cam.f, l = cam.l,r = cam.r, b = cam.b, t = cam.t, nx = cam.sizeX, ny = cam.sizeY;

        double Mper[4][4] = {{(2*n) / (r-l), 0, (r+l) / (r-l), 0},
                            {0, (2*n) / (t-b), (t+b) / (t-b), 0},
                            {0, 0, -((f+n) / (f-n)), -((2*f*n) / (f-n))},
                            {0, 0, -1, 0}};

        double Mvp[4][4] = {{nx/2.0, 0, 0, (nx-1) / 2.0 },
                           {0, ny/2.0, 0, (ny-1) / 2.0 },
                           {0, 0, 0.5, 0.5 },
                           {0, 0, 0, 1}};

         multiplyMatrixWithMatrix(Mcam, R, E);

        double depth[numberOfVertices];

        for(int i = 1; i <= numberOfVertices; i++) {
                double r[4];
                Vec3 cur_vertex = changedVertices[i];
                double cur_vertex_vec[4] = {cur_vertex.x, cur_vertex.y, cur_vertex.z, 1};

                //Camera transformation
                multiplyMatrixWithVec4d(r, Mcam, cur_vertex_vec);
                double cam_transformed_vertex[4] = {r[0], r[1], r[2], r[3]};

                //Perspective transformation
                multiplyMatrixWithVec4d(r, Mper, cam_transformed_vertex );
                changedVertices[i] = {r[0], r[1], r[2]};
                depth[i] = r[3];
        }

        for (int a = 1; a <= numberOfVertices; a++) {
                changedVertices2[a] = changedVertices[a];
        }

        Model changedModels[numberOfModels];
        for(int i = 0; i < numberOfModels; i++) {
                changedModels[i] = models[i];
        }

        //Backface culling
        if (backfaceCullingSetting) {
                for(int i = 0; i < numberOfModels; i++) {
                        Model model = changedModels[i];

                        for(int j = 0; j < model.numberOfTriangles; j++) {
                                Triangle cur_tri = (model.triangles)[j];
                                Vec3 a = changedVertices[cur_tri.vertexIds[0]];
                                Vec3 b = changedVertices[cur_tri.vertexIds[1]];
                                Vec3 c = changedVertices[cur_tri.vertexIds[2]];
                                Vec3 middle_point = {(a.x + b.x + c.x)/3,
                                                     (a.y + b.y + c.y)/3,
                                                     (a.z + b.z + c.z)/3};

                                Vec3 view_vec = subtractVec3(middle_point, {0,0,0});
                                Vec3 normal_vec = crossProductVec3(subtractVec3(c,b), subtractVec3(a,b));
                                if(dotProductVec3(view_vec, normal_vec) <= 0) {
                                        changedModels[i].triangles[j].vertexIds[0] = -99;

                                }
                        }
                }
        }

        for(int i = 1; i <= numberOfVertices; i++) {
                //Perspective division
                double cam_transformed_vertex2[4] = {changedVertices2[i].x / depth[i],changedVertices2[i].y / depth[i], changedVertices2[i].z / depth[i], 1};
                double r[4];
                //Viewport transformation
                multiplyMatrixWithVec4d(r, Mvp, cam_transformed_vertex2 );
                changedVertices2[i].x = r[0];
                changedVertices2[i].y = r[1];
                changedVertices2[i].z = r[2];
        }

      //Rasterization
      for(int i = 0; i < numberOfModels; i++) {
              Model model = changedModels[i];
              for(int j = 0; j < model.numberOfTriangles; j++) {
                      Triangle cur_tri = (model.triangles)[j];
                      //cur_vertex_ids = triangle's vertices
                      int *cur_vertex_ids = cur_tri.vertexIds;

                      double x0 = changedVertices2[cur_vertex_ids[0]].x;
                      double x1 = changedVertices2[cur_vertex_ids[1]].x;
                      double x2 = changedVertices2[cur_vertex_ids[2]].x;

                      double y0 = changedVertices2[cur_vertex_ids[0]].y;
                      double y1 = changedVertices2[cur_vertex_ids[1]].y;
                      double y2 = changedVertices2[cur_vertex_ids[2]].y;

                      if (cur_vertex_ids[0] != -99) {
                              Color color0 = colors[cur_vertex_ids[0]];
                              Color color1 = colors[cur_vertex_ids[1]];
                              Color color2 = colors[cur_vertex_ids[2]];

                              Vec3 color_vec0 = {color0.r, color0.g, color0.b};
                              Vec3 color_vec1 = {color1.r, color1.g, color1.b};
                              Vec3 color_vec2 = {color2.r, color2.g, color2.b};

                              if (model.type == 0) {
                                      if (x1 > x0)
                                              midpointLine( x1, x0, y1, y0, color_vec0, color_vec1, cam );
                                      else
                                              midpointLine( x0, x1, y0, y1, color_vec1, color_vec0, cam );

                                      if (x2 > x0)
                                              midpointLine( x2, x0, y2, y0, color_vec0, color_vec2, cam );
                                      else
                                              midpointLine( x0, x2, y0, y2, color_vec2, color_vec0, cam );

                                      if (x2 > x1)
                                              midpointLine( x2, x1, y2, y1, color_vec1, color_vec2, cam );
                                      else
                                              midpointLine( x1, x2, y1, y2, color_vec2, color_vec1, cam );
                              }

                              else {
                                      double xmin = min3(x0 ,x1, x2);
                                      double ymin = min3(y0 ,y1, y2);
                                      double xmax = max3(x0 ,x1, x2);
                                      double ymax = max3(y0 ,y1, y2);

                                      for(int k = xmin; k <= xmax; k++) {
                                              for(int l = ymin; l <= ymax; l++) {
                                                      double alpha = fxy(k, l, x1, x2, y1, y2) / fxy(x0, y0, x1, x2, y1, y2);
                                                      double beta  = fxy(k, l, x2, x0, y2, y0) / fxy(x1, y1, x2, x0, y2, y0);
                                                      double gama  = fxy(k, l, x0, x1, y0, y1) / fxy(x2, y2, x0, x1, y0, y1);
                                                      if(alpha >= 0 && beta >= 0 && gama >= 0) {

                                                              Vec3 color_m0 = multiplyVec3WithScalar(color_vec0,alpha);
                                                              Vec3 color_m1 = multiplyVec3WithScalar(color_vec1,beta);
                                                              Vec3 color_m2 = multiplyVec3WithScalar(color_vec2,gama);

                                                              Vec3 new_color = addVec3(addVec3(color_m0, color_m1), color_m2);
                                                              if (k < nx && l < ny && k > 0 && l > 0) {
                                                                      image[k][l].r = (int)new_color.x;
                                                                      image[k][l].g = (int)new_color.y;
                                                                      image[k][l].b = (int)new_color.z;
                                                              }
                                                      }
                                              }
                                      }
                              }
                      }
              }
      }
}


int main(int argc, char **argv) {
    int i, j;

    if (argc < 2) {
        std::cout << "Usage: ./rasterizer <scene file> <camera file>" << std::endl;
        return 1;
    }

    // read camera and scene files
    readSceneFile(argv[1]);
    readCameraFile(argv[2]);

    image = 0;

    for (i = 0; i < numberOfCameras; i++) {

        // allocate memory for image
        if (image) {
			for (j = 0; j < cameras[i].sizeX; j++) {
		        delete image[j];
		    }

			delete[] image;
		}

        image = new Color*[cameras[i].sizeX];

        if (image == NULL) {
            std::cout << "ERROR: Cannot allocate memory for image." << std::endl;
            exit(1);
        }

        for (j = 0; j < cameras[i].sizeX; j++) {
            image[j] = new Color[cameras[i].sizeY];
            if (image[j] == NULL) {
                std::cout << "ERROR: Cannot allocate memory for image." << std::endl;
                exit(1);
            }
        }


        // initialize image with basic values
        initializeImage(cameras[i]);

        // do forward rendering pipeline operations
        forwardRenderingPipeline(cameras[i]);

        // generate PPM file
        writeImageToPPMFile(cameras[i]);

        // Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
        // Notice that os_type is not given as 1 (Ubuntu) or 2 (Windows), below call doesn't do conversion.
        // Change os_type to 1 or 2, after being sure that you have ImageMagick installed.
        convertPPMToPNG(cameras[i].outputFileName, 99);
    }

    return 0;

}
