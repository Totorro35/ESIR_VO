/**
 * @file av-v0.cpp
 * @author Thomas LEMETAYER (thomas.lemetayer.35@gmail.com)
 * @author Corentin SALAUN (corentin.salaun@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2019-10-15
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#define VP_TRACE

//! \example tutorial-viewer.cpp
//! [Include display]
#include <visp3/gui/vpDisplayD3D.h>
#include <visp3/gui/vpDisplayGDI.h>
#include <visp3/gui/vpDisplayGTK.h>
#include <visp3/gui/vpDisplayX.h>
#include <visp3/gui/vpDisplayOpenCV.h>
#include <visp3/core/vpPoint.h>
#include <visp3/core/vpPlane.h>
#include <visp3/gui/vpPlot.h>
#include <visp3/core/vpMeterPixelConversion.h>
#include <visp3/core/vpCameraParameters.h>

#include <visp3/core/vpExponentialMap.h>
#include <visp3/core/vpVelocityTwistMatrix.h>



//! [Include display]
//! [Include io]
#include <visp3/io/vpImageIo.h>
//! [Include io]

using namespace std ;


void display(vpCameraParameters& cam, vpImage<unsigned char> &I, vpColVector &x, vpColVector &xd)
{
    for (int i = 0 ; i < x.getRows()/2 ; i++)
    {
        vpImagePoint u,ud ;
        vpMeterPixelConversion::convertPoint(cam,x[2*i],x[2*i+1],u) ;
        vpDisplay::displayPoint(I,u,vpColor::red, 2) ;

        vpMeterPixelConversion::convertPoint(cam,xd[2*i],xd[2*i+1],ud) ;
        vpDisplay::displayCircle(I,ud,10,vpColor::green) ;

    }

    vpDisplay::flush(I) ;
}


// Projection d'un point 3D sur le plane image  X(3), x(2)
void project(vpColVector &X, vpColVector &x)
{

}

// Changement de repere bX(3), aX(3), aTb est une matrice homogène
/**
 * @brief 
 * 
 * @param bX vecteur exprime dans le repere a
 * @param aTb Matrice de transformation du repere a au repere b
 * @param aX Output vector X exprime dans le repere a
 */
void changeFrame(const vpColVector &bX, const vpHomogeneousMatrix &aTb, vpColVector &aX)
{

}



// Calcul de la matrice d'interaction d'un point 2D
void computeInteractionMatrix(vpColVector &cX, double x, double y, vpMatrix &Lx)
{


}




void tp2DVisualServoingOnePoint()
{

    //-------------------------------------------------------------
    // Mise en oeuvre des courbes



    vpPlot plot(4, 700, 700, 100, 200, "Curves...");


    char title[40];
    strncpy( title, "||e||", 40 );
    plot.setTitle(0,title);
    plot.initGraph(0,1);

    strncpy( title, "x-xd", 40 );
    plot.setTitle(1, title);
    plot.initGraph(1,2);

    strncpy( title, "camera velocity", 40 );
    plot.setTitle(2, title);
    plot.initGraph(2,6);


    strncpy( title, "Camera position", 40 );
    plot.setTitle(3, title);
    plot.initGraph(3,6);

    //-------------------------------------------------------------
    // Affichage des images
    vpImage<unsigned char> I(400,600) ;
    vpDisplayX d ;
    d.init(I) ;
    vpDisplay::display(I);
    vpCameraParameters cam(400,400,300,200) ;

    //-------------------------------------------------------------

    //Definition de la scene

    //positions initiale (à tester)

    //Definition de la scene
    vpHomogeneousMatrix cTw (0,0,1,  0,0,0) ;

    //position of the point in the world frame
    vpColVector wX(3) ;
    wX[0] = 0.5 ; wX[1] = 0.2 ; wX[2] = -0.5 ;


    // a chaque fois que vous verez size metter la bonne taille à la place
    int size ;

    vpColVector e(size) ; //


    // position courante, position desiree
    vpColVector x(size), xd(size) ;
    //matrice d'interaction
    vpMatrix Lx(size,size) ;

    // position desirée  (au centre de l'image x=0, y=0)


    // vitesse de la camera
    vpColVector v(size) ;
    double lambda = 0.1 ;
    int iter = 0 ;
    while (fabs(e.sumSquare()) > 1e-6)
    {
        // calcul de la position des points dans l'image en fonction de cTw
       // instancier x


        //calcul de l'erreur

        // Calcul de la matrice d'interaction


        //calcul de la loi de commande v= ...

        v = 0 ;


        // Ne pas modifier la suite
        //mise a jour de la position de la camera
        cTw = vpExponentialMap::direct(v).inverse()* cTw ;

        cout << "iter "<< iter <<" : "<< e.t() << endl ;

        iter++ ;



        //mise a jour des courbes
        vpPoseVector ctw(cTw) ;
        plot.plot(0,0,iter, e.sumSquare()) ;
        plot.plot(1,iter, e) ;
        plot.plot(2,iter, v) ;
        plot.plot(3,iter, ctw) ;
        //mise a jour de l'image
        display(cam,I,x,xd) ;


    }

    // sauvegarde des courbes

    plot.saveData(0,"e.txt","#");
    plot.saveData(1,"error.txt","#");
    plot.saveData(2,"v.txt","#");
    plot.saveData(3,"p.txt","#");

    int a ; cin >> a ;
    // sauvegarde de l'image finale
    {
        vpImage<vpRGBa>  Irgb ;
        vpDisplay::getImage(I,Irgb) ;
        vpImageIo::write(Irgb,"1pt.jpg") ;
    }
    cout << "Clicker sur l'image pour terminer" << endl ;
    vpDisplay::getClick(I) ;
}

void tp2DVisualServoingFourPoint()
{



    //-------------------------------------------------------------
    // Mise en oeuvre des courbes

    vpPlot plot(4, 700, 700, 100, 200, "Curves...");


    char title[40];
    strncpy( title, "||e||", 40 );
    plot.setTitle(0,title);
    plot.initGraph(0,1);

    strncpy( title, "x-xd", 40 );
    plot.setTitle(1, title);
    plot.initGraph(1,8);

    strncpy( title, "camera velocity", 40 );
    plot.setTitle(2, title);
    plot.initGraph(2,6);


    strncpy( title, "camera position", 40 );
    plot.setTitle(3, title);
    plot.initGraph(3,6);

    //-------------------------------------------------------------
    // Affichage des images
    vpImage<unsigned char> I(400,600) ;
    vpDisplayX d ;
    d.init(I) ;
    vpDisplay::display(I);
    vpCameraParameters cam(400,400,300,200) ;

    //-------------------------------------------------------------


    //positions initiale (à tester)
    vpHomogeneousMatrix cTw (-0.2, -0.1, 1.3,
                             vpMath::rad(10), vpMath::rad(20), vpMath::rad(30) ) ;
    //  vpHomogeneousMatrix cTw (0.2,0.1,1.3,  0,0,vpMath::rad(5)) ;
    //   vpHomogeneousMatrix cTw (0,0,1,  0,0,vpMath::rad(45)) ;
    //  vpHomogeneousMatrix cTw (0, 0, 1,  0, 0, vpMath::rad(90)) ;
    //   vpHomogeneousMatrix cTw (0, 0, 1,  0, 0, vpMath::rad(180)) ;

    // position finale
    vpHomogeneousMatrix cdTw (0,0,1,  0,0,0) ;



    // position des point dans le repere monde Rw
    vpColVector wX[4] ;
    for (int i = 0 ; i < 4 ; i++) wX[i].resize(3) ;

    double M = 0.3 ;
    wX[0][0] = -M     ; wX[0][1] = -M     ; wX[0][2] = 0 ;
    wX[1][0] = M ;      wX[1][1] = -M     ; wX[1][2] = 0 ;
    wX[2][0] = M ;      wX[2][1] =  M;      wX[2][2] = 0 ;
    wX[3][0] =  -M    ; wX[3][1] =  M;      wX[3][2] = 0 ;

    int size ;
    vpColVector e(size) ; //

    vpColVector x(size), xd(size) ;

    //initialisation de la position désire des points dans l'image en fonction de cdTw


    vpColVector v(size) ;
    double lambda = 0.1 ;
    int iter = 0 ;


    while (fabs(e.sumSquare()) > 1e-16)
    {

        // calcul de la position des points dans l'image en fonction de cTw


        // Calcul de la matrice d'interaction
        vpMatrix Lx(size,size) ;


        //calcul de l'erreur



        //calcul de la loi de commande


        //mise a jour de la position de la camera
        cTw = vpExponentialMap::direct(v).inverse()* cTw ;

        cout << "iter "<< iter <<" : "<< e.t() << endl ;
        iter++ ;

       //mise a jour des courbes
        vpPoseVector ctw(cTw) ;
        plot.plot(0,0,iter, e.sumSquare()) ;
        plot.plot(1,iter, e) ;
        plot.plot(2,iter, v) ;
        plot.plot(3,iter, ctw) ;
        //mise a jour de l'image
        display(cam,I,x,xd) ;

    }

    // sauvegarde des courbes
    plot.saveData(0,"e.txt","#");
    plot.saveData(1,"error.txt","#");
    plot.saveData(2,"v.txt","#");
    plot.saveData(3,"p.txt","#");

    // sauvegarde de l'image finale
    {
        vpImage<vpRGBa>  Irgb ;
        vpDisplay::getImage(I,Irgb) ;
        vpImageIo::write(Irgb,"4pt.jpg") ;
    }
    cout << "Clicker sur l'image pour terminer" << endl ;
    vpDisplay::getClick(I) ;
}

/*
void
computeError3D(...)
{


}

void
computeInteractionMatrix3D(...)
{


}
*/

void tp3DVisualServoing()
{


    vpTRACE("begin" ) ;

    vpPlot plot(4, 700, 700, 100, 200, "Curves...");


    char title[40];
    strncpy( title, "||e||", 40 );
    plot.setTitle(0,title);
    plot.initGraph(0,1);

    strncpy( title, "x-xd", 40 );
    plot.setTitle(1, title);
    plot.initGraph(1,6);

    strncpy( title, "camera velocity", 40 );
    plot.setTitle(2, title);
    plot.initGraph(2,6);


    strncpy( title, "Camera position", 40 );
    plot.setTitle(3, title);
    plot.initGraph(3,6);



    //Definition de la scene
    vpHomogeneousMatrix cTw (-0.2, -0.1, 1.3,
                             vpMath::rad(10), vpMath::rad(20), vpMath::rad(30) ) ;
    vpHomogeneousMatrix cdTw (0,0,1,  0,0,0) ;



    int size ;
    vpColVector e(size) ; //

    vpMatrix Lx(size,size) ;

    vpColVector v(size) ;
    double lambda = 0.1 ;
    int iter = 0 ;
    while (fabs(e.sumSquare()) > 1e-6)
    {

        // Calcul de l'erreur
        //computeError3D(...) ;

        // Calcul de la matrice d'interaction
        //  computeInteractionMatrix3D(...) ;

        //        Calcul de la loi de commande


        // Mis à jour de la position de la camera
        cTw = vpExponentialMap::direct(v).inverse()* cTw ;

        cout << "iter "<< iter <<" : "<< e.t() << endl ;

        iter++ ;

    //mis a jour de courbes
        vpPoseVector crw(cTw) ;
        plot.plot(0,0,iter, e.sumSquare()) ;
        plot.plot(1,iter, e) ;
        plot.plot(2,iter, v) ;
        plot.plot(3,iter,crw) ;


    }

// sauvegarde des courbes
    plot.saveData(0,"e.txt","#");
    plot.saveData(1,"error.txt","#");
    plot.saveData(2,"v.txt","#");
    plot.saveData(3,"p.txt","#");

    int a ; cin >> a ;

}


void tp2DVisualServoingFourPointMvt()
{

/*
    reprendre votre code pour 4 points
    au debut de la boucle d'av vous considerer que les points sont en mouvement (o
            while (fabs(e.sumSquare()) > 1e-16)
            {
                // les points sont animés d'un mouvement de 1cm/s en x dans Rw
                for (int i = 0 ; i < 4 ; i++) wX[i][0] += 0.01 ;
            }

*/
}




int main(int argc, char** argv)
{

  //  tp2DVisualServoingOnePoint() ;
    // tp2DVisualServoingFourPoint() ;
    //tp3DVisualServoing() ;
      tp2DVisualServoingFourPointMvt() ;

}
