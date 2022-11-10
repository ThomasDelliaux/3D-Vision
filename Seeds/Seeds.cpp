// Imagine++ project
// Project:  Seeds
// Author:   Pascal Monasse
// Student: TODO

#include <Imagine/Images.h>
#include <queue>
#include <string>
#include <iostream>
using namespace Imagine;
using namespace std;

// Default data
const char *DEF_im1=srcPath("face0.jpg"), *DEF_im2=srcPath("face1.jpg");
//static int dmin=-30, dmax=-7; // Min and max disparities
static int dmin=10, dmax=55;
/// Min NCC for a seed
static const float nccSeed=0.95f;

/// Radius of patch for correlation
static const int win=(9-1)/2;
/// To avoid division by 0 for constant patch
static const float EPS=0.1f;

/// A seed
struct Seed {
    Seed(int x0, int y0, int d0, float ncc0)
    : x(x0), y(y0), d(d0), ncc(ncc0) {}
    int x,y, d;
    float ncc;
};

/// Order by NCC
bool operator<(const Seed& s1, const Seed& s2) {
    return (s1.ncc<s2.ncc);
}

/// 4-neighbors
static const int dx[]={+1,  0, -1,  0};
static const int dy[]={ 0, -1,  0, +1};

/// Display disparity map
static Image<Color> displayDisp(const Image<int>& disp, Window W, int subW) {
    Image<Color> im(disp.width(), disp.height());
    for(int j=0; j<disp.height(); j++)
        for(int i=0; i<disp.width(); i++) {
            if(disp(i,j)<dmin || disp(i,j)>dmax)
                im(i,j) = CYAN;
            else {
                int g = 255*(disp(i,j)-dmin)/(dmax-dmin);
                im(i,j)= Color(g,g,g);
            }
        }
    setActiveWindow(W,subW);
    display(im);
    showWindow(W,subW);
    return im;
}

/// Show 3D window
static void show3D(const Image<Color>& im, const Image<int>& disp) {
#ifdef IMAGINE_OPENGL // Imagine++ must have been built with OpenGL support...
    // Intrinsic parameters given by Middlebury website
    const float f=3740;
    const float d0=-200; // Doll images cropped by this amount
    const float zoom=2; // Half-size images, should double measured disparity
    const float B=0.160; // Baseline in m
    FMatrix<float,3,3> K(0.0f);
    K(0,0)=-f/zoom; K(0,2)=disp.width()/2;
    K(1,1)= f/zoom; K(1,2)=disp.height()/2;
    K(2,2)=1.0f;
    K = inverse(K);
    K /= K(2,2);
    std::vector<FloatPoint3> pts;
    std::vector<Color> col;
    for(int j=0; j<disp.height(); j++)
        for(int i=0; i<disp.width(); i++)
            if(dmin<=disp(i,j) && disp(i,j)<=dmax) {
                float z = B*f/(zoom*disp(i,j)+d0);
                FloatPoint3 pt((float)i,(float)j,1.0f);
                pts.push_back(K*pt*z);
                col.push_back(im(i,j));
            }
    Mesh mesh(&pts[0], pts.size(), 0,0,0,0,VERTEX_COLOR);
    mesh.setColors(VERTEX, &col[0]);
    Window W = openWindow3D(512,512,"3D");
    setActiveWindow(W);
    showMesh(mesh);
#else
    std::cout << "No 3D: Imagine++ not built with OpenGL support" << std::endl;
#endif
}

/// Correlation between patches centered on (i1,j1) and (i2,j2). The values
/// m1 or m2 are subtracted from each pixel value.
static float correl(const Image<byte>& im1, int i1,int j1,float m1,
                    const Image<byte>& im2, int i2,int j2,float m2) {
    float dist=0.0f;
    //Those variables are used for store the different sums, required for the computation of the NCC.
    float upper_sum=0.0f;
    float lower_sum_i1 = 0.0f;
    float lower_sum_i2 = 0.0f;

    //We pass through each pixel of the patches.
    for(int k{-win};k<=win;++k){
        for(int l{-win};l<=win;++l){
            upper_sum+=(im1(i1+k,j1+l)-m1)*(im2(i2+k,j2+l)-m2);
            lower_sum_i1+=pow(im1(i1+k,j1+l)-m1,2);
            lower_sum_i2+=pow(im2(i2+k,j2+l)-m2,2);

        }
    }
    //add EPS for avoir division by 0.
    dist = upper_sum/(sqrt(lower_sum_i1)*sqrt(lower_sum_i2)+EPS);
    return dist;
}

/// Sum of pixel values in patch centered on (i,j).
static float sum(const Image<byte>& im, int i, int j) {
    float s=0.0f;
    //Pass through each pixel of the patch.
    for(int k{-win};k<=win;++k){
        for(int l{-win};l<=win;++l){
            s+=(im(i+k,j+l));

        }
    }
    return s;
}

/// Centered correlation of patches of size 2*win+1.
static float ccorrel(const Image<byte>& im1,int i1,int j1,
                     const Image<byte>& im2,int i2,int j2) {
    float m1 = sum(im1,i1,j1);
    float m2 = sum(im2,i2,j2);
    int w = 2*win+1;
    return correl(im1,i1,j1,m1/(w*w), im2,i2,j2,m2/(w*w));
}

/// Compute disparity map from im1 to im2, but only at points where NCC is
/// above nccSeed. Set to true the seeds and put them in Q.
static void find_seeds(Image<byte> im1, Image<byte> im2,
                       float nccSeed,
                       Image<int>& disp, Image<bool>& seeds,
                       std::priority_queue<Seed>& Q) {
    disp.fill(dmin-1);
    seeds.fill(false);
    while(! Q.empty())
        Q.pop();

    const int maxy = std::min(im1.height(),im2.height());
    const int refreshStep = (maxy-2*win)*5/100;
    for(int y=win; y+win<maxy; y++) {
        if((y-win-1)/refreshStep != (y-win)/refreshStep)
            std::cout << "Seeds: " << 5*(y-win)/refreshStep <<"%\r"<<std::flush;

        for(int x=win; x+win<im1.width(); x++) {
            //Store the best ncc found for a patch. bestd is the disparity associated to bestNcc.
            float bestNcc = -5.0f;
            int bestd = 0;

            //Test each disparity between dmin and dmax included for find the ncc of the patch.
            for(int d = dmin;d<=dmax;++d){
                //Check if the moved patch belongs to im2.
                if(x+d-win>=0 && x+d+win<im2.width()){

                    float ncc = ccorrel(im1,x,y,im2,x+d,y);
                    //If we find a better ncc we replace bestNcc and bestd.
                    if(ncc>=bestNcc){
                        bestNcc = ncc;
                        bestd = d;
                    }


                 }
             }
            //Check if the ncc found for a patch is superior to a threshold for the seeds algorithm.
            if(bestNcc>=nccSeed){

                disp(x,y)=bestd;
                seeds(x,y) = true;
                Seed newSeed(x,y,bestd,bestNcc);
                Q.push(newSeed);

            }
        }
    }
    std::cout << std::endl;
}

/// Propagate seeds
static void propagate(Image<byte> im1, Image<byte> im2,
                      Image<int>& disp, Image<bool>& seeds,
                      std::priority_queue<Seed>& Q) {
    const int maxy = std::min(im1.height(),im2.height());
    while(! Q.empty()) {
        Seed s=Q.top();
        Q.pop();
        for(int i=0; i<4; i++) {
            int x=s.x+dx[i], y=s.y+dy[i];
            if(0<=x-win && x+win<im1.width() && 0<=y-win && y+win<maxy &&
               ! seeds(x,y)) {

                float ncc = -1.0f, nccPlus1 = -1.0f, nccMinus1 = -1.0f;
                //Check if the patch still belongs to the image 2.
                if(0<=x+s.d-win && x+s.d+win<im2.width()){
                    ncc = ccorrel(im1,x,y,im2,x+s.d,y);
                }
                //Check if the patch still belongs to the image 2 and if d+1 is in [dmin,dmax]
                if(0<=x+s.d+1-win && x+s.d+1+win<im1.width() && s.d<dmax){
                     nccPlus1 = ccorrel(im1,x,y,im2,x+s.d+1,y);
                }
                //Check if the patch still belongs to the image 2 and if d-1 is in [dmin,dmax]
                if(0<=x+s.d-1-win && x+s.d-1+win<im1.width() && s.d>dmin){
                     nccMinus1 = ccorrel(im1,x,y,im2,x+s.d-1,y);
                }
                //Find the disparity which maximise the ncc between d,d+1,d-1. And then update the seeds map, the disparity map,
                //and put the current point in the seeds array.
                if(ncc>=nccPlus1 && ncc>=nccMinus1){
                    seeds(x,y)=true;
                    disp(x,y)=s.d;
                    Seed newSeed(x,y,s.d,ncc);
                    Q.push(newSeed);
                }
                else if(nccPlus1>=ncc && nccPlus1>=nccMinus1){
                    seeds(x,y)=true;
                    disp(x,y)=s.d+1;
                    Seed newSeed(x,y,s.d+1,nccPlus1);
                    Q.push(newSeed);

                }
                else{
                    seeds(x,y)=true;
                    disp(x,y)=s.d-1;
                    Seed newSeed(x,y,s.d-1,nccMinus1);
                    Q.push(newSeed);

                }

            }
        }
    }
}

int main(int argc, char* argv[]) {
    if(argc!=1 && argc!=5) {
        cerr << "Usage: " << argv[0] << " im1 im2 dmin dmax" << endl;
        return 1;
    }
    const char *im1=DEF_im1, *im2=DEF_im2;
    if(argc>1) {
        im1 = argv[1]; im2=argv[2]; dmin=stoi(argv[3]); dmax=stoi(argv[4]);
    }
    // Load and display images
    Image<Color> I1, I2;
    if(!load(I1,im1) || !load(I2,im2)) {
        cerr<< "Error loading image files" << endl;
        return 1;
    }
    //I1=I1.getSubImage(0,0,400,400);
    //I2=I2.getSubImage(0,0,400,400);
    std::string names[5]={"image 1","image 2","dense","seeds","propagation"};
    Window W = openComplexWindow(I1.width(), I1.height(), "Seeds propagation",
                                 5, names);
    setActiveWindow(W,0);
    display(I1,0,0);
    setActiveWindow(W,1);
    display(I2,0,0);

    Image<int> disp(I1.width(), I1.height());
    Image<bool> seeds(I1.width(), I1.height());
    std::priority_queue<Seed> Q;

    // Dense disparity
    find_seeds(I1, I2, -1.0f, disp, seeds, Q);
    save(displayDisp(disp,W,2), srcPath("0dense.png"));

    // Only seeds
    find_seeds(I1, I2, nccSeed, disp, seeds, Q);
    save(displayDisp(disp,W,3), srcPath("1seeds.png"));

    // Propagation of seeds
    propagate(I1, I2, disp, seeds, Q);
    save(displayDisp(disp,W,4), srcPath("2final.png"));

    // Show 3D (use shift click to animate)
    show3D(I1,disp);

    endGraphics();
    return 0;
}
