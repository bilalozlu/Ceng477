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

/*P3
# 0
0 255
255

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
double min3(double n1, double n2, double n3){
        if (n2>=n3 && n1>=n3 ){
                return n3;
        }
        else if (n1>=n2 && n3>=n2){
                return n2;
        }
        else{
                return n1;
        }
}

double max3(double n1, double n2, double n3){
        if (n1>=n2 && n1>=n3 ){
                return n1;
        }
        else if (n2>=n1 && n2>=n3){
                return n2;
        }
        else{
                return n3;
        }
}

double fxy(double x, double y, double x0, double x1, double y0, double y1){
        double res = x*(y0-y1) + y*(x1-x0) + x0*y1 -y0*x1;
        // std::cout<<"res: "<<res<<std::endl;
        return res;
}

void translation(const Translation &S, const Vec3 &in, Vec3 &out ) {

        double r[4];
        double m1[4][4] = {{1,0,0,S.tx},
                        {0,1,0,S.ty},
                        {0,0,1,S.tz},
                        {0,0,0,1}};
        double m2[4] = {in.x, in.y, in.z,1};

        multiplyMatrixWithVec4d(r,m1,m2);

        out = {r[0],r[1],r[2]};

}

void scaling(const Scaling &S , const Vec3 &in , Vec3 &out) {

        double r[4];
        double m1[4][4] = {{S.sx,0,0,0},
                        {0,S.sy,0,0},
                        {0,0,S.sz,0},
                        {0,0,0,1}};
        double m2[4] = {in.x,in.y,in.z,1};

        multiplyMatrixWithVec4d(r,m1,m2);

        out = {r[0],r[1],r[2]};

}

void rotation_transform(const Rotation &S, double out[4][4] ) {


        double r1[4][4];
        Vec3 V = {-S.uy,S.ux,0};
        Vec3 W = crossProductVec3({S.ux,S.uy,S.uz},V);
        double m1[4][4] = {{S.ux,S.uy,S.uz,0},
                        {V.x,V.y,V.z,0},
                        {W.x,W.y,W.z,0},
                        {0,0,0,1}};

        double m1_inverse[4][4] = {{S.ux,V.x,W.x,0},
                        {S.uy,V.y,W.y,0},
                        {S.uz,V.z,W.z,0},
                        {0,0,0,1}};

        double rota1[4][4] = {{1,0,0,0},
                        {0,cos(S.angle * PI / 180.0),-sin(S.angle * PI / 180.0),0},
                        {0,sin(S.angle * PI / 180.0),cos(S.angle * PI / 180.0),0},
                        {0,0,0,1}};

        multiplyMatrixWithMatrix(r1,m1_inverse,rota1);
        multiplyMatrixWithMatrix(out,r1,m1);
        // std::cout<<"111111111111111"<<std::endl;
        // for(int i=0;i<4;i++) {
        //     for(int j=0;j<4;j++) {
        //         std::cout<<out[i][j]<<" ";
        //     }
        //     std::cout<<"\n";
        // }
        // std::cout<<"\n\n\n\n";
}

void forwardRenderingPipeline(Camera cam) {

    std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<< all the vertices: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;
    for(int i=1;i<=numberOfVertices;i++) {
        std::cout<<vertices[i].x<<" "<<vertices[i].y<<" "<<vertices[i].z<<std::endl;
    }
    std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<< all the vertices: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;

    std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<< all the modelIds, modelTypes, : <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;
    std::cout<<"number of models: "<<numberOfModels<<std::endl;
    for(int i=0;i<numberOfModels;i++) {
        std::cout<<models[i].modelId<<" "<<models[i].type<<std::endl;
    }
    std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<< all the modelIds, modelTypes, : <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;

        bool isTransformed[1000];
        Vec3 changedVertices[numberOfVertices+1];
        Vec3 changedVertices2[numberOfVertices+1];

        for (int a = 1; a<= numberOfVertices; a++){
                changedVertices[a] = vertices[a];
        }

        //Modeling transformations
        for (int j = 0; j < numberOfModels; j++) {
                Model model = models[j];
                for(int i = 0; i < model.numberOfTransformations; i++) {
                        for(int k = 0; k < model.numberOfTriangles; k++) {
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
                                                        // std::cout<<"changedVertices before: "<<changedVertices[cur_vertex_id].x<<" "
                                                        // <<changedVertices[cur_vertex_id].y<<" "<<changedVertices[cur_vertex_id].z<<std::endl;
                                                        translation(translations[trans_id], vertex, changedVertices[cur_vertex_id]);
                                                        // std::cout<<"changedVertices after: "<<changedVertices[cur_vertex_id].x<<" "
                                                        // <<changedVertices[cur_vertex_id].y<<" "<<changedVertices[cur_vertex_id].z<<std::endl;
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
                                        double res[4];
                                        rotation_transform(rotations[trans_id], theMatrix);
                                        for (int l = 0; l < 3; l++) {
                                                int cur_vertex_id = cur_tri.vertexIds[l];
                                                if (!isTransformed[cur_vertex_id]) {
                                                        Vec3 vertex = changedVertices[cur_vertex_id];
                                                        isTransformed[cur_vertex_id] = true;
                                                        double theVector[4] = {vertex.x,vertex.y,vertex.z,1.0};
                                                        multiplyMatrixWithVec4d(res, theMatrix, theVector);
                                                        changedVertices[cur_vertex_id] = {res[0],res[1],res[2]};
                                                }
                                        }
                                }
                        }
                        else {
                                std::cout << "No such transformation";
                                return;
                        }
                }
        }

        double Mcam[4][4];
        //  = {{cam.u.x,cam.u.y, cam.u.z,-(cam.u.x*cam.pos.x+cam.u.y*cam.pos.y+cam.u.z*cam.pos.z)},
        //                      {cam.v.x, cam.v.y, cam.v.z,-(cam.v.x*cam.pos.x+cam.v.y*cam.pos.y+cam.v.z*cam.pos.z)},
        //                      {cam.w.x, cam.w.y, cam.w.z,-(cam.w.x*cam.pos.x+cam.w.y*cam.pos.y+cam.w.z*cam.pos.z)},
        //                      {0, 0, 0, 1.0}};

        double R[4][4] = {{cam.u.x, cam.u.y, cam.u.z, 0},
                        {cam.v.x, cam.v.y, cam.v.z, 0},
                        {cam.w.x, cam.w.y, cam.w.z, 0},
                        {0, 0, 0, 1}};

        double E[4][4] = {{1, 0, 0, -cam.pos.x},
                        {0, 1, 0, -cam.pos.y},
                        {0, 0, 1, -cam.pos.z},
                        {0, 0, 0, 1}};

        double n = cam.n, f = cam.f, l = cam.l,r = cam.r, b = cam.b, t = cam.t, nx = cam.sizeX, ny = cam.sizeY;

        // std::cout<<"Camera Vars-------------------------------------------------------------------------------------"<<std::endl;
        // std::cout<<"n: "<<n<<" "<<" f: "<<f<<" l: "<<l<<" r: "<<
        // r<<" b: "<<b<<" t: "<<t<<"nx: "<<nx<<"ny: "<<ny<<std::endl;
        // std::cout<<"Camera Vars-------------------------------------------------------------------------------------"<<std::endl;

        double Mper[4][4] = {{(2*n)/(r-l), 0, (r+l)/(r-l), 0},
                          {0, (2*n)/(t-b), (t+b)/(t-b), 0},
                          {0, 0, -((f+n)/(f-n)), -((2*f*n)/(f-n))},
                          {0, 0, -1, 0}};

        double Mvp[4][4] = {{nx/2.0, 0, 0, (nx-1)/2.0 },
                          {0, ny/2.0, 0, (ny-1)/2.0 },
                          {0, 0, 0.5, 0.5 },
                          {0, 0, 0, 1}};

         multiplyMatrixWithMatrix(Mcam, R, E);


        for(int i=1;i<=numberOfVertices;i++) {
            double r[4];
            Vec3 cur_vertex = changedVertices[i];
            double cur_vertex_vec[4] = {cur_vertex.x,cur_vertex.y,cur_vertex.z,1};
            // std::cout<<"before camera: "<<cur_vertex_vec[0]<<"      "<< cur_vertex_vec[1]<<"     "
            // <<"      "<<cur_vertex_vec[2]<<"       "<<cur_vertex_vec[3]<<std::endl;
            //Camera transformation
            multiplyMatrixWithVec4d(r, Mcam, cur_vertex_vec);

            double cam_transformed_vertex[4] = {r[0], r[1], r[2], r[3]};
            // std::cout<<"after camera: "<<r[0]<<"      "<< r[1]<<"     "<<"      "<<r[2]<<"       "<<r[3]<<std::endl;

            // perspective transformation
            multiplyMatrixWithVec4d(r, Mper, cam_transformed_vertex );
            //Projection division
            // std::cout << " after perspective transformation: "<< r[0]<<
            // "             "<<r[1]<<"             "<<r[2]<<"      "<<r[3] <<std::endl;

            double cam_transformed_vertex2[4] = {r[0]/r[3], r[1]/r[3], r[2]/r[3], 1};

            //Viewport transformation
            // std::cout << " after project division: "<< cam_transformed_vertex2[0]<<"      "<<cam_transformed_vertex2[1]<<"      "<<cam_transformed_vertex2[2]<< std::endl;

            multiplyMatrixWithVec4d(r, Mvp, cam_transformed_vertex2 );
            changedVertices[i] = {r[0], r[1], r[2]};

            std::cout << " after viewport: "<< r[0]<<"      "<<r[1]<<"      "<<r[2]<< std::endl;

        }


        for(int i=1; i<=numberOfVertices; i++) {
            changedVertices2[i] = changedVertices[i];
        }

        //Backface culling
      //   for(int i=0;i<numberOfModels;i++) {
      //           Model model = models[i];
      //           for(int j=0;j<model.numberOfTriangles;j++) {
      //                   Triangle cur_tri = (model.triangles)[j];
      //                   Vec3 a = changedVertices[cur_tri.vertexIds[0]];
      //                   Vec3 b = changedVertices[cur_tri.vertexIds[1]];
      //                   Vec3 c = changedVertices[cur_tri.vertexIds[2]];
      //                   Vec3 middle_point = {(a.x + b.x + c.x)/3,
      //                                   (a.y + b.y + c.y)/3,
      //                                   (a.z + b.z + c.z)/3};
      //
      //           Vec3 view_vector = subtractVec3(middle_point, cam.pos);
      //           Vec3 normal_vector = crossProductVec3(subtractVec3(c,b), subtractVec3(a,b));
      //           Vec3 minus = {-1, -1, -1};
      //
      //           changedVertices2[cur_tri.vertexIds[0]] = dotProductVec3(view_vector, normal_vector) > 0 ? changedVertices[cur_tri.vertexIds[0]] : minus;
      //           changedVertices2[cur_tri.vertexIds[1]] = dotProductVec3(view_vector, normal_vector) > 0 ? changedVertices[cur_tri.vertexIds[1]] : minus;
      //           changedVertices2[cur_tri.vertexIds[2]] = dotProductVec3(view_vector, normal_vector) > 0 ? changedVertices[cur_tri.vertexIds[2]] : minus;
      //   }
      // }

//       for ( int i = 0; i < numberOfVertices; i++ ){
//               double r[4];
//               //Projection division
//               double cam_transformed_vertex2[4] = {changedVertices[i].x/w[i], changedVertices[i].y/w[i], changedVertices[i].z/w[i], 1};
//               //Viewport transformation
//               //multiplyMatrixWithVec4d(r, Mvp, cam_transformed_vertex2 );
//               changedVertices[i] = {cam_transformed_vertex2[0], cam_transformed_vertex2[1], cam_transformed_vertex2[2]};
// }

      // for ( int i = 0; i < numberOfVertices; i++ ){
      //         std::cout << "c" << changedVertices[i].x << " " << changedVertices[i].y << " " << changedVertices[i].z  << std::endl;
      //         std::cout << "o" << vertices[i].x << " " << vertices[i].y << " " << vertices[i].z  << std::endl;
      // }

      // for(int i = 0;i<numberOfModels;i++) {
      //     Model model = models[i];
      //     for(int j=0;j<model.numberOfTriangles;j++) {
      //             Triangle cur_tri = (model.triangles)[j];
      //             int *cur_vertex_ids = cur_tri.vertexIds;
      //             colors[changedVertices[cur_vertex_ids[i+1]].colorId] = {10,10,10};
      //     }
      // }

      std::cout<<"colors"<<std::endl;
      for(int i=1;i<=numberOfColors;i++) {

          std::cout<<colors[i].r<<" "<<colors[i].g<<" "<<colors[i].r<<std::endl;
      }


      //Rasterization
      for(int i = 0;i<numberOfModels;i++) {
              Model model = models[i];
              // std::cout<<"number of Triangles: "<<model.numberOfTriangles<<std::endl;
              for(int j=0;j<model.numberOfTriangles;j++) {

                      Triangle cur_tri = (model.triangles)[j];
                      // cur_vertex_ids means the triangle's vertices
                      int *cur_vertex_ids = cur_tri.vertexIds;

                      std::cout<<cur_vertex_ids[0]<<" "<<cur_vertex_ids[1]<<" "<<cur_vertex_ids[2]<<std::endl;

                      int x0 = changedVertices[cur_vertex_ids[0]].x;
                      int x1 = changedVertices[cur_vertex_ids[1]].x;
                      int x2 = changedVertices[cur_vertex_ids[2]].x;

                      int y0 = changedVertices[cur_vertex_ids[0]].y;
                      int y1 = changedVertices[cur_vertex_ids[1]].y;
                      int y2 = changedVertices[cur_vertex_ids[2]].y;
                      //noktasal renkler. Sonradan sil
                      // (image[x0][y0]) = colors[cur_vertex_ids[0]];
                      // (image[x1][y1]) = colors[cur_vertex_ids[1]];
                      // (image[x2][y2]) = colors[cur_vertex_ids[2]];

                      //noktasal renklerin hepsi siyah. Sonradan sil
                      (image[x0][y0]).r = 0; (image[x0][y0]).g = 0; (image[x0][y0]).b = 0;
                      (image[x1][y1]).r = 0; (image[x1][y1]).g = 0; (image[x1][y1]).b = 0;
                      (image[x2][y2]).r = 0; (image[x2][y2]).g = 0; (image[x2][y2]).b = 0;


                      // std::cout<<"----------------------------- Two dimensional triangle values -------------------------"<<std::endl;
                      // std::cout<<x0<<" "<<y0<<std::endl;
                      // std::cout<<x1<<" "<<y1<<std::endl;
                      // std::cout<<x2<<" "<<y2<<std::endl;
                      // std::cout<<"----------------------------- Two dimensional triangle values -------------------------"<<std::endl;

                      // int xmin = min3(x0 ,x1,x2);
                      // int ymin = min3(y0 ,y1,y2);
                      // int xmax = max3(x0 ,x1, x2);
                      // int ymax = max3(y0 ,y1, y2);

                      // for(int k=xmin;k<=xmax;k++) {
                      //         for(int l=ymin;l<=ymax;l++) {
                      //
                      //                 double alpha = -1;
                      //                 double fVal = fxy(x0,y0,x1,x2,y1,y2);
                      //                 if(fVal >= 0.000001){
                      //                     alpha = fxy(k,l,x1,x2,y1,y2)/fVal;
                      //                 }
                      //
                      //                 double beta = -1;
                      //                 double fVal2 = fxy(x1,y1,x1,x2,y1,y2);
                      //                 if(fVal2 >= 0.000001){
                      //                     beta = fxy(k,l,x2,x0,y2,y0)/fVal2;
                      //                 }
                      //
                      //                 double gama = -1;
                      //                 double fVal3 = fxy(x2,y2,x1,x2,y1,y2);
                      //                 if(fVal3 >= 0.000001){
                      //                     gama = fxy(k,l,x0,x1,y0,y1)/fVal3;
                      //                 }

                                      // double x, double y, double x0, double x1, double y0, double y1

                                      // double beta  = fxy(k,l,x2,x0,y2,y0)/fxy(x1,y1,x1,x2,y1,y2);
                                      // double gama  = fxy(k,l,x0,x1,y0,y1)/fxy(x2,y2,x1,x2,y1,y2);
                                      //
                                      // if(fVal == 0 || fVal2 == 0 || fVal3 == 0 || (alpha >=0 && beta >=0 && gama >=0)) {
                                      //     //
                                      //     // std::cout<<"--------------------- Barycentric coordinates ------------------"<<std::endl;
                                      //     // std::cout<<"alpha: "<<alpha<<" beta: "<<beta<<" gama: "<<gama<<std::endl;
                                      //     // std::cout<<"---------------- Barycentric coordinates ---------------"<<std::endl;
                                      //
                                      //         Color color0 = colors[cur_vertex_ids[0]];
                                      //         Color color1 = colors[cur_vertex_ids[1]];
                                      //         Color color2 = colors[cur_vertex_ids[2]];
                                      //         std::cout<<changedVertices[cur_vertex_ids[0]].colorId<<std::endl;
                                      //         // if(changedVertices[cur_vertex_ids[0]].colorId !=0)
                                      //         //   std::cout<<"ids: "<<changedVertices[cur_vertex_ids[0]].colorId<<" "<<changedVertices[cur_vertex_ids[1]].colorId<<" "<<changedVertices[cur_vertex_ids[2]].colorId<<std::endl;
                                      //
                                      //
                                      //           // std::cout<<"ids22: "<<cur_vertex_ids[0]<<" "<<cur_vertex_ids[1]<<" "<<cur_vertex_ids[2]<<std::endl;
                                      //
                                      //       std::cout<<"ids: "<<color1.r<<" " <<color1.g<<" "<<color1.b<<std::endl;
                                      //       // std::cout<<"ids: "<<color2.r<<" " <<color2.g<<" "<<color2.b<<std::endl;
                                      //
                                      //         Vec3 color_vec0 = {color0.r, color0.g, color0.b};
                                      //         Vec3 color_vec1 = {color1.r, color1.g, color1.b};
                                      //         Vec3 color_vec2 = {color2.r, color2.g, color2.b};
                                      //
                                      //         Vec3 color_m0 = multiplyVec3WithScalar(color_vec0,alpha);
                                      //         Vec3 color_m1 = multiplyVec3WithScalar(color_vec1,beta);
                                      //         Vec3 color_m2 = multiplyVec3WithScalar(color_vec2,gama);
                                      //
                                      //         Vec3 new_color = addVec3(addVec3(color_m0, color_m1), color_m2);
                                      //         // std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<       new_color      >>>>>>>>>>>>>>>>>>>>>"<<std::endl;
                                      //         // std::cout<<new_color.x<<" "<<new_color.y<<" "<<new_color.z<<std::endl;
                                      //         // std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<       new_color      >>>>>>>>>>>>>>>>>>>>>"<<std::endl;
                                      //         image[k][l].r = new_color.x;
                                      //         image[k][l].g = new_color.y;
                                      //         image[k][l].b = new_color.z;
                                      //
                                      // }
                      //         }
                      // }
              }
      }

      for(int i = 1;i<=numberOfVertices; i++){
            vertices[i] = changedVertices2[i];
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

    for (i = 0; i < 1; i++) {

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
