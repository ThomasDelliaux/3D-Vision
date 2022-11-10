// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure

struct Match {
    float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color,2> I1, Image<Color,2> I2,
              vector<Match>& matches) {
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0,0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(),0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0*100.0;
    for(size_t i=0; i < feats1.size(); i++) {
        SIFTDetector::Feature f1=feats1[i];
        for(size_t j=0; j < feats2.size(); j++) {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if(d < MAX_DISTANCE) {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}

//Compute a fundamental matrix from eight matching points.
FMatrix<float,3,3> computeRandF(vector<Match> & eightPoints){


    FMatrix<float,9,9> A;
    A.fill(0);


    //Normalization matrix for normalized the coordinates of our points. In order to make the svd in a efficient way.
    FMatrix<float,3,3> normMatrice;
    normMatrice.fill(0);
    normMatrice(0,0) = 1e-3;
    normMatrice(1,1) = 1e-3;
    normMatrice(2,2) = 1;


    //Fill the matrice A from the course formula.
    for(int i{0}; i<8;++i){
        FVector<float,3> x_t;
        FVector<float,3> x_t_p;


        x_t[0] = eightPoints[i].x1, x_t[1] = eightPoints[i].y1, x_t[2] = 1;
        x_t_p[0] = eightPoints[i].x2, x_t_p[1] = eightPoints[i].y2, x_t_p[2] = 1;

        //x_t correspond to the x tilde in the cours dans x_t_p correspond to the x tilde prime in the course.
        x_t = normMatrice*x_t;
        x_t_p = normMatrice*x_t_p;


        A(i,0) = x_t[0]*x_t_p[0];
        A(i,1) = x_t[0]*x_t_p[1];
        A(i,2) = x_t[0];
        A(i,3) = x_t[1]*x_t_p[0];
        A(i,4) = x_t[1]*x_t_p[1];
        A(i,5) = x_t[1];
        A(i,6) = x_t_p[0];
        A(i,7) = x_t_p[1];
        A(i,8) = 1;
    }

    //Perfom the svd on A.
    FMatrix<float,9,9> U, Vt;
    FVector<float,9> S;
    FVector<float,9> f;
    FMatrix<float,3,3> F;
    F.fill(0);

    svd(A,U,S,Vt);

    f=transpose(Vt).getCol(8);

    F(0,0) = f[0], F(0,1) = f[1], F(0,2) = f[2];
    F(1,0) = f[3], F(1,1) = f[4], F(1,2) = f[5];
    F(2,0) = f[6], F(2,1) = f[7], F(2,2) = f[8];


    //Ensure that the rank of the matrix F is 2.
    FMatrix<float,3,3> U_f, V_f;
    FVector<float,3> S_f;

    svd(F,U_f,S_f,V_f);

    S_f[2] = 0;

    F = U_f*Diagonal(S_f)*V_f;
    F = normMatrice*F*normMatrice;

    return F;


}

//Does the same thing as computeRandF but this time we use the set of inliers found by RANSAC algorithm instead of the eight random points.
FMatrix<float,3,3> refineWithLeastSquarred(vector<int> & bestInliers,vector<Match> & matches){
    int NBestInliers = (int)bestInliers.size();

    Matrix<float> A(NBestInliers,9);
    A.fill(0.);


    //Normalization matrix for normalized the coordinates of our points. In order to make the svd in a efficient way.
    FMatrix<float,3,3> normMatrice;
    normMatrice.fill(0);
    normMatrice(0,0) = 1e-3;
    normMatrice(1,1) = 1e-3;
    normMatrice(2,2) = 1;


    //Fill the matrice A from the course formula.
    for(int i{0}; i<NBestInliers;++i){
        int indexBestInliers = bestInliers[i];
        FVector<float,3> x_t;
        FVector<float,3> x_t_p;


        x_t[0] = matches[indexBestInliers].x1, x_t[1] = matches[indexBestInliers].y1, x_t[2] = 1;
        x_t_p[0] = matches[indexBestInliers].x2, x_t_p[1] = matches[indexBestInliers].y2, x_t_p[2] = 1;

        //x_t correspond to the x tilde in the cours dans x_t_p correspond to the x tilde prime in the course.
        x_t = normMatrice*x_t;
        x_t_p = normMatrice*x_t_p;


        A(i,0) = x_t[0]*x_t_p[0];
        A(i,1) = x_t[0]*x_t_p[1];
        A(i,2) = x_t[0];
        A(i,3) = x_t[1]*x_t_p[0];
        A(i,4) = x_t[1]*x_t_p[1];
        A(i,5) = x_t[1];
        A(i,6) = x_t_p[0];
        A(i,7) = x_t_p[1];
        A(i,8) = 1;
    }

    //Perfom the svd on A.
    Matrix<float> U(NBestInliers,NBestInliers), Vt(9,9);
    Vector<float> S(NBestInliers);
    Vector<float> f(9);
    FMatrix<float,3,3> F;
    F.fill(0);

    svd(A,U,S,Vt);

    f=Vt.getRow(8);

    F(0,0) = f[0], F(0,1) = f[1], F(0,2) = f[2];
    F(1,0) = f[3], F(1,1) = f[4], F(1,2) = f[5];
    F(2,0) = f[6], F(2,1) = f[7], F(2,2) = f[8];


    //Ensure that the rank of the matrix F is 2.
    FMatrix<float,3,3> U_f, V_f;
    FVector<float,3> S_f;

    svd(F,U_f,S_f,V_f);

    S_f[2] = 0;

    F = U_f*Diagonal(S_f)*V_f;
    F = normMatrice*F*normMatrice;

    return F;
}
// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3> computeF(vector<Match>& matches) {
    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    int Niter=100000; // Adjusted dynamically
    FMatrix<float,3,3> bestF;

    //It contains the indexes of the largest set of inliers found by RANSAC.
    vector<int> bestInliers;

    int n = matches.size();
    vector<Match> eightPoints(8);

    int k{0};
    while(k<=Niter){

        //Pick randomly 8 matching points in the matches.

        for(int j{0}; j<=7; ++j){

            eightPoints[j] = matches[rand()%(int)n];
        }

        FMatrix<float,3,3> F;
        F = computeRandF(eightPoints);

        //Perfom the RANSAC algorithm.

        //Use this vector to store the index of the detected inliers.
        vector<int> currentInliers;

        for(int i{0};i<n; ++i){
            FVector<float,3> x(3);

            //Correspond to the x prime in the course.
            FVector<float,3> x_p(3);

            x[0] = matches[i].x1, x[1] = matches[i].y1, x[2] = 1;
            x_p[0] = matches[i].x2, x_p[1] = matches[i].y2, x_p[2]=1;

            //Epipolar line associated to x.
            FVector<float,3> Ft_x = transpose(F)*x;

            float distToEpipolar = abs((Ft_x*x_p))/sqrt(pow(Ft_x[0],2)+pow(Ft_x[1],2));

            if (distToEpipolar < distMax){
                currentInliers.push_back(i);
            }
        }

        if (currentInliers.size()>bestInliers.size()){
            bestInliers=currentInliers;

            //Check that there are at least 10% of the number of matches in bestInliers before updating the number of iterations, to avoid numerical issue.
            if((bestInliers.size())>n*0.1){
                Niter = ceil((log(BETA)/log(1-pow(((float)bestInliers.size()/(float)n),8))));

            }
        }


        ++k;





    }

    //Use the best intliers found during the RANSAC algorithm for refine our matrix.
    bestF = refineWithLeastSquarred(bestInliers, matches);

    cout<<bestInliers.size()<<" inliers detected !"<<endl;
    // Updating matches with inliers only
    vector<Match> all=matches;
    matches.clear();
    for(size_t i=0; i<bestInliers.size(); i++)
        matches.push_back(all[bestInliers[i]]);
    return bestF;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float,3,3>& F) {

    while(true) {
        int x,y;
        if(getMouse(x,y) == 3)
            break;

        FVector<float,3> point;
        FVector<float,3> epipolarLine;

        int w1{I1.width()};
        int w2{I2.width()};

        point[0] = x, point[1] = y, point[2] = 1;

        //Check if (x,y) belongs to the image 1 or the image 2.
        if(x<I1.width()){
            drawCircle(x,y,5,BLUE);
            point[0] = x, point[1] = y, point[2] = 1;
            epipolarLine = transpose(F)*point;
            drawLine(w1,-epipolarLine[2]/epipolarLine[1],w1+w2,-(epipolarLine[2]+w2*epipolarLine[0])/epipolarLine[1],RED);
        }

        if(x>=I1.width()){
            drawCircle(x,y,5,BLUE);

            //x-w1 it's to put (0,0) in the upper left corner of the image 2.
            point[0] = x-w1, point[1] = y, point[2] = 1;
            epipolarLine = F*point;
            drawLine(0,-epipolarLine[2]/epipolarLine[1],w1,-(epipolarLine[2]+w1*epipolarLine[0])/epipolarLine[1],RED);
        }


    }
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    const char* s1 = argc>1? argv[1]: srcPath("im1.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("im2.jpg");

    // Load and display images
    Image<Color,2> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2*w, I1.height());
    display(I1,0,0);
    display(I2,w,0);

    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    const int n = (int)matches.size();
    cout << " matches: " << n << endl;
    drawString(100,20,std::to_string(n)+ " matches",RED);

    cout<<"Click on the image to compute the fundamental matrix, and display the largest set of inliers found by RANSAC algorithm."<<endl;
    click();
    
    FMatrix<float,3,3> F = computeF(matches);
    cout << "F="<< endl << F;

    // Redisplay with matches
    display(I1,0,0);
    display(I2,w,0);
    for(size_t i=0; i<matches.size(); i++) {
        Color c(rand()%256,rand()%256,rand()%256);
        fillCircle(matches[i].x1+0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2+w, matches[i].y2, 2, c);        
    }
    drawString(100, 20, to_string(matches.size())+"/"+to_string(n)+" inliers", RED);
    cout<<"Click on the image to continue."<<endl;
    click();

    cout<<"Click on a point to display the associated epipolar line. Use the right click for close the image."<<endl;
    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}
