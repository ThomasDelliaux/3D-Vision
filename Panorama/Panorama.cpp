// Imagine++ project
// Project:  Panorama
// Author:   Pascal Monasse
// Date:     2013/10/08

#include <Imagine/Graphics.h>
#include <Imagine/Images.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <sstream>
using namespace Imagine;
using namespace std;

// Record clicks in two images, until right button click
void getClicks(Window w1, Window w2,
               vector<IntPoint2>& pts1, vector<IntPoint2>& pts2) {

    //Display number of each image for practical purpose.
    setActiveWindow(w1);
    drawString(10,15,"Image 1",RED);

    setActiveWindow(w2);
    drawString(10,15,"Image 2",RED);


    int whichClick;
    IntPoint2 clickedPoint;
    int pointNumber{0};

    setActiveWindow(w1);
    showWindow(w1);


    cout<<"Please select at least 4 points on the image 1, and then press the right click."<<endl;

    //Request to click on points as long as the right click is not used or the number of selected points is less than 4.
    //This ensures that the user has selected at lest 4 points.

    while((whichClick = getMouse(clickedPoint))!=3 || pointNumber<=3){
        if(whichClick==1){
        pts1.push_back(clickedPoint);
        fillCircle(clickedPoint,2,RED);
        ++pointNumber;
        cout<<"You have selected "<<pointNumber<<" point."<<endl;
        }
        else{
            cout<<"You have to select at least "<<4<<" points"<<endl;
        }

        }

    //After that the program switch to the second image and asks the same thing.
    pointNumber = 0;

    setActiveWindow(w2);
    showWindow(w2);

    cout<<"Please select at least 4 points on the image 2, and then press the right click."<<endl;


    while((whichClick = getMouse(clickedPoint))!=3 || pointNumber<=3){

        if(whichClick){
        pts2.push_back(clickedPoint);
        fillCircle(clickedPoint,2,RED);
        ++pointNumber;
        cout<<"You have selected "<<pointNumber<<" point."<<endl;
        }
        else{
            cout<<"You have to select at least "<<4<<" points"<<endl;
        }
    }
}

// Return homography compatible with point matches
Matrix<float> getHomography(const vector<IntPoint2>& pts1,
                            const vector<IntPoint2>& pts2) {
    size_t n = min(pts1.size(), pts2.size());
    if(n<4) {
        cout << "Not enough correspondences: " << n << endl;
        return Matrix<float>::Identity(3);
    }

    //Fill the matrix with 0s because I observed that it's increase the stability of my program.

    Matrix<double> A(2*n,8);
    A.fill(0);
    Vector<double> B(2*n);
    B.fill(0);

    //Fill the matrix A and B according to the course formula.
    for(int i{0};i<int(n);++i){
        int x_i_1= pts1[i].x();
        int y_i_1= pts1[i].y();
        int x_i_2= pts2[i].x();
        int y_i_2= pts2[i].y();

        A(2*i,0)=x_i_1;
        A(2*i,1)=y_i_1;
        A(2*i,2) = 1;
        A(2*i,6)=-x_i_2*x_i_1;
        A(2*i,7) = -x_i_2*y_i_1;

        A(2*i+1,3) = x_i_1;
        A(2*i+1,4) = y_i_1;
        A(2*i+1,5) = 1;
        A(2*i+1,6) = -y_i_2*x_i_1;
        A(2*i+1,7) =-y_i_1*y_i_2;

        B[2*i] = x_i_2;
        B[2*i+1] = y_i_2;
   }


    //Solve the linear system to recover the homography's coefficents.
    B = linSolve(A, B);
    Matrix<float> H(3, 3);
    H(0,0)=B[0]; H(0,1)=B[1]; H(0,2)=B[2];
    H(1,0)=B[3]; H(1,1)=B[4]; H(1,2)=B[5];
    H(2,0)=B[6]; H(2,1)=B[7]; H(2,2)=1;

    // Sanity check
    for(size_t i=0; i<n; i++) {
        float v1[]={(float)pts1[i].x(), (float)pts1[i].y(), 1.0f};
        float v2[]={(float)pts2[i].x(), (float)pts2[i].y(), 1.0f};
        Vector<float> x1(v1,3);
        Vector<float> x2(v2,3);
        x1 = H*x1;
        cout << x1[1]*x2[2]-x1[2]*x2[1] << ' '
             << x1[2]*x2[0]-x1[0]*x2[2] << ' '
             << x1[0]*x2[1]-x1[1]*x2[0] << endl;
    }
    return H;
}

// Grow rectangle of corners (x0,y0) and (x1,y1) to include (x,y)
void
growTo(float& x0, float& y0, float& x1, float& y1, float x, float y) {
    if(x<x0) x0=x;
    if(x>x1) x1=x;
    if(y<y0) y0=y;
    if(y>y1) y1=y;    
}

// Panorama construction
void panorama(const Image<Color,2>& I1, const Image<Color,2>& I2,
              Matrix<float> H) {

    Vector<float> v(3);
    float x0=0, y0=0, x1=I2.width(), y1=I2.height();

    v[0]=0; v[1]=0; v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=I1.width(); v[1]=0; v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=I1.width(); v[1]=I1.height(); v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=0; v[1]=I1.height(); v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    cout << "x0 x1 y0 y1=" << x0 << ' ' << x1 << ' ' << y0 << ' ' << y1<<endl;

    Image<Color> I(int(x1-x0), int(y1-y0));
    setActiveWindow( openWindow(I.width(), I.height()) );
    I.fill(WHITE);

    //Use the pull back method to perfom the panorama algorithm.
    Matrix<float> inverseH{inverse(H)};

    //I distinguish the coordinates of the plane of image 2 and image 1.
    Vector<float> coordPlane2(3);
    Vector<float> coordPlane1(3);


    //Go through each pixel and check if this pixel belongs to the image 2 and if its image by the inverted homography belongs to image 1.
    for(int i{0};i<I.width();++i){
        for(int j{0};j<I.height();++j){

            //Allow to stock the information about if the pixel belongs to image 1 or 2.
            bool inI1=false;
            bool inI2=false;

            //Allows to start at the "true" first pixel on the plane 2(the plane of the image 2).
            coordPlane2[0]=i+x0, coordPlane2[1]=j+y0, coordPlane2[2]=1;

            //Use the inverted homography to project the pixel on the plane 1(the plane of the image 1).
            coordPlane1 = inverseH*coordPlane2; coordPlane1/=coordPlane1[2];



            //Check if the pixel(in plane 1 coordinates) belongs to the image 1.
            if((0<=coordPlane1[0] && coordPlane1[0]<I1.width()) && (coordPlane1[1]>=0 && coordPlane1[1]<I1.height())){
                inI1=true;
            }

            //Check if the pixel(in plane 2 coordinates) belongs to the image 2.
            if((0<=coordPlane2[0] && coordPlane2[0]<I2.width()) && (coordPlane2[1]>=0 && coordPlane2[1]<I2.height())){
                I(i,j)=I2(coordPlane2[0],coordPlane2[1]);
                inI2=true;

            }

            //If the pixel belongs to image 1 and image 2, make a mean of the colors.
            if(inI1 && inI2){
                RGB<int> deux (2,2,2);
               I(i,j)= I1(coordPlane1[0],coordPlane1[1])/(unsigned char)2;
               I(i,j)+=I2(coordPlane2[0],coordPlane2[1])/(unsigned char)2;
            }

            //If the pixels belongs to the image 1/image 2, then fill the pixel with its color in the image 1/image 2.
            else if(inI1){
                I(i,j)=I1(coordPlane1[0],coordPlane1[1]);
            }
            else if(inI2){
                I(i,j)=I2(coordPlane2[0],coordPlane2[1]);
            }




       }
    }
    display(I,0,0);
}

// Main function
int main(int argc, char* argv[]) {
    const char* s1 = argc>1? argv[1]: srcPath("image0006.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("image0007.jpg");

    // Load and display images
    Image<Color> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load the images" << endl;
        return 1;
    }
    Window w1 = openWindow(I1.width(), I1.height(), s1);
    display(I1,0,0);
    Window w2 = openWindow(I2.width(), I2.height(), s2);
    setActiveWindow(w2);
    display(I2,0,0);

    // Get user's clicks in images
    vector<IntPoint2> pts1, pts2;
    getClicks(w1, w2, pts1, pts2);

    vector<IntPoint2>::const_iterator it;
    cout << "pts1="<<endl;
    for(it=pts1.begin(); it != pts1.end(); it++)
        cout << *it << endl;
    cout << "pts2="<<endl;
    for(it=pts2.begin(); it != pts2.end(); it++)
        cout << *it << endl;

    // Compute homography
    Matrix<float> H = getHomography(pts1, pts2);
    cout << "H=" << H/H(2,2);

    // Apply homography
    panorama(I1, I2, H);

    endGraphics();
    return 0;
}
