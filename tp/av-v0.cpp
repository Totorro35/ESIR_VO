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

#define CUDAVITE

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

#include <visp3/io/vpImageIo.h>

using namespace std;

/**
 * @brief Display an the View
 * 
 * @param cam 
 * @param I 
 * @param x 
 * @param xd 
 */
void display(vpCameraParameters &cam, vpImage<unsigned char> &I, vpColVector &x, vpColVector &xd)
{
    for (int i = 0; i < x.getRows() / 2; i++)
    {
        vpImagePoint u, ud;
        vpMeterPixelConversion::convertPoint(cam, x[2 * i], x[2 * i + 1], u);
        vpDisplay::displayPoint(I, u, vpColor::red, 2);

        vpMeterPixelConversion::convertPoint(cam, xd[2 * i], xd[2 * i + 1], ud);
        vpDisplay::displayCircle(I, ud, 10, vpColor::green);
    }

    vpDisplay::flush(I);
}

/**
 * @brief Projection d'un point 3D sur le plane image  X(3), x(2)
 * 
 * @param X 
 * @param x 
 */
void project(vpColVector &X, vpColVector &x)
{
    x.resize(2);
    x[0] = X[0] / X[2];
    x[1] = X[1] / X[2];
}

/**
 * @brief Changement de repere bX(3), aX(3), aTb est une matrice homogène
 * 
 * @param bX Vecteur dans le repère B
 * @param aTb Matrice de transformation
 * @param aX Output vecteur dans le repère A
 */
void changeFrame(const vpColVector &bX, const vpHomogeneousMatrix &aTb, vpColVector &aX)
{
    aX.resize(bX.size());
    vpColVector bXHomogene(bX.size() + 1), aXHomogene(bX.size() + 1);
    for (size_t i = 0; i < bXHomogene.size() - 1; i++)
    {
        bXHomogene[i] = bX[i];
    }
    bXHomogene[aXHomogene.size() - 1] = 1.0;

    aXHomogene = aTb * bXHomogene;

    for (size_t i = 0; i < aXHomogene.size() - 1; i++)
    {
        aX[i] = aXHomogene[i] / aXHomogene[aXHomogene.size() - 1];
    }
}

/**
 * @brief Calcul de la matrice d'interaction d'un point 2D
 * 
 * @param cX 
 * @param x 
 * @param y 
 * @param Lx 
 */
void computeInteractionMatrix(vpColVector &cX, double x, double y, vpMatrix &Lx)
{
    double invZ = 1 / cX[2];
    Lx.resize(2, 6);
    Lx[0][0] = -invZ;
    Lx[0][1] = 0;
    Lx[0][2] = x * invZ;
    Lx[0][3] = x * y;
    Lx[0][4] = -(1 + x * x);
    Lx[0][5] = y;

    Lx[1][0] = 0;
    Lx[1][1] = -invZ;
    Lx[1][2] = y * invZ;
    Lx[1][3] = 1 + y * y;
    Lx[1][4] = -x * y;
    Lx[1][5] = -x;
}

/**
 * @brief Executable for 3.2
 *  Asservissement visuel sur 1 point 2D
 * 
 */
void tp2DVisualServoingOnePoint()
{

    //-------------------------------------------------------------
    // Mise en oeuvre des courbes

    vpPlot plot(4, 700, 700, 100, 200, "Curves...");

    char title[40];
    strncpy(title, "||e||", 40);
    plot.setTitle(0, title);
    plot.initGraph(0, 1);

    strncpy(title, "x-xd", 40);
    plot.setTitle(1, title);
    plot.initGraph(1, 2);

    strncpy(title, "camera velocity", 40);
    plot.setTitle(2, title);
    plot.initGraph(2, 6);

    strncpy(title, "Camera position", 40);
    plot.setTitle(3, title);
    plot.initGraph(3, 6);

    //-------------------------------------------------------------
    // Affichage des images
    vpImage<unsigned char> I(400, 600);
    vpDisplayX d;
    d.init(I);
    vpDisplay::display(I);
    vpCameraParameters cam(400, 400, 300, 200);

    //-------------------------------------------------------------

    //Definition de la scene

    //positions initiale (à tester)

    //Definition de la scene
    vpHomogeneousMatrix cTw(0, 0, 1, 0, 0, 0);

    //position of the point in the world frame
    vpColVector wX(3);
    wX[0] = 0.5;
    wX[1] = 0.2;
    wX[2] = -0.5;

    vpColVector e(2); //
    e = 1;

    // position courante, position desiree
    vpColVector x(2), xd(2);
    //matrice d'interaction
    vpMatrix Lx(2, 6);

    // position desirée  (au centre de l'image x=0, y=0)
    xd[0] = 0;
    xd[1] = 0;

    // vitesse de la camera
    vpColVector v(6);
    double lambda = 0.1;
    int iter = 0;
    while (fabs(e.sumSquare()) > 1e-6)
    {
        // calcul de la position des points dans l'image en fonction de cTw
        // instancier x

        vpColVector cX;
        changeFrame(wX, cTw, cX);
        project(cX, x);

        //calcul de l'erreur
        e = x - xd;

        // Calcul de la matrice d'interaction
        computeInteractionMatrix(cX, x[0], x[1], Lx);

        //calcul de la loi de commande v= ...

        v = -lambda * Lx.pseudoInverse() * e;

        // Ne pas modifier la suite
        //mise a jour de la position de la camera
        cTw = vpExponentialMap::direct(v).inverse() * cTw;

        cout << "iter " << iter << " : " << e.t() << endl;

        iter++;

        //mise a jour des courbes
        vpPoseVector ctw(cTw);
        plot.plot(0, 0, iter, e.sumSquare());
        plot.plot(1, iter, e);
        plot.plot(2, iter, v);
        plot.plot(3, iter, ctw);
        //mise a jour de l'image
        display(cam, I, x, xd);
    }

    // sauvegarde des courbes

    plot.saveData(0, "e.txt", "#");
    plot.saveData(1, "error.txt", "#");
    plot.saveData(2, "v.txt", "#");
    plot.saveData(3, "p.txt", "#");

    int a;
    cin >> a;
    // sauvegarde de l'image finale
    {
        vpImage<vpRGBa> Irgb;
        vpDisplay::getImage(I, Irgb);
        vpImageIo::write(Irgb, "1pt.jpg");
    }
    cout << "Clicker sur l'image pour terminer" << endl;
    vpDisplay::getClick(I);
}

/**
 * @brief Executable for 3.3
 *  Asservissement visuel sur 4 points 2D
 * 
 */
void tp2DVisualServoingFourPoint()
{

    //-------------------------------------------------------------
    // Mise en oeuvre des courbes

    vpPlot plot(4, 700, 700, 100, 200, "Curves...");

    char title[40];
    strncpy(title, "||e||", 40);
    plot.setTitle(0, title);
    plot.initGraph(0, 1);

    strncpy(title, "x-xd", 40);
    plot.setTitle(1, title);
    plot.initGraph(1, 8);

    strncpy(title, "camera velocity", 40);
    plot.setTitle(2, title);
    plot.initGraph(2, 6);

    strncpy(title, "camera position", 40);
    plot.setTitle(3, title);
    plot.initGraph(3, 6);

    //-------------------------------------------------------------
    // Affichage des images
    vpImage<unsigned char> I(400, 600);
    vpDisplayX d;
    d.init(I);
    vpDisplay::display(I);
    vpCameraParameters cam(400, 400, 300, 200);

    //-------------------------------------------------------------

    //positions initiale (à tester)
    //vpHomogeneousMatrix cTw(0, 0, 1.3, 0, 0, 0);
    //vpHomogeneousMatrix cTw(-0.2, -0.1, 1.3,vpMath::rad(10), vpMath::rad(20), vpMath::rad(30));
    //vpHomogeneousMatrix cTw(0.2, 0.1, 1.3, 0, 0, vpMath::rad(5));
    //vpHomogeneousMatrix cTw(0, 0, 1, 0, 0, vpMath::rad(45));
    //vpHomogeneousMatrix cTw(0, 0, 1, 0, 0, vpMath::rad(90));
    vpHomogeneousMatrix cTw(0, 0, 1, 0, 0, vpMath::rad(180));

    // position finale
    vpHomogeneousMatrix cdTw(0, 0, 1, 0, 0, 0);

    // position des point dans le repere monde Rw
    vpColVector wX[4];
    for (int i = 0; i < 4; i++)
        wX[i].resize(3);

    double M = 0.3;
    wX[0][0] = -M;
    wX[0][1] = -M;
    wX[0][2] = 0;
    wX[1][0] = M;
    wX[1][1] = -M;
    wX[1][2] = 0;
    wX[2][0] = M;
    wX[2][1] = M;
    wX[2][2] = 0;
    wX[3][0] = -M;
    wX[3][1] = M;
    wX[3][2] = 0;

    vpColVector e(8); //
    e = 1;

    vpColVector x(8), xd(8);

    //initialisation de la position désire des points dans l'image en fonction de cdTw

    for (size_t i = 0; i < 4; i++)
    {
        vpColVector cdX;
        changeFrame(wX[i], cdTw, cdX);
        vpColVector _xd;
        project(cdX, _xd);
        xd[i * 2 + 0] = _xd[0];
        xd[i * 2 + 1] = _xd[1];
    }

    vpColVector v(6);
    double lambda = 0.1;
    int iter = 0;
    while (fabs(e.sumSquare()) > 1e-16)
    {
        std::cout << "----------------------------" << std::endl;
        std::cout << "Itération : " << iter << std::endl;
        vpMatrix Lx(8, 6);

        // calcul de la position des points dans l'image en fonction de cTw
        for (size_t i = 0; i < 4; i++)
        {
            vpColVector _cX;
            changeFrame(wX[i], cTw, _cX);

            vpColVector _x;
            project(_cX, _x);

            x[i * 2 + 0] = _x[0];
            x[i * 2 + 1] = _x[1];

            vpColVector _cdX;
            changeFrame(wX[i], cdTw, _cdX);
            vpColVector _xd;
            project(_cdX, _xd);

            vpMatrix _Lx;
            computeInteractionMatrix(_cX, _x[0], _x[1], _Lx);
            //computeInteractionMatrix(_cdX, _xd[0], _xd[1], _Lx);

            for (size_t j = 0; j < 6; j++)
            {
                Lx[i * 2 + 0][j] = _Lx[0][j];
                Lx[i * 2 + 1][j] = _Lx[1][j];
            }
        }

        //calcul de l'erreur
        e = x - xd;
        //calcul de la loi de commande
        v = -lambda * Lx.pseudoInverse() * e;

        //mise a jour de la position de la camera
        cTw = vpExponentialMap::direct(v).inverse() * cTw;

        cout << "iter " << iter << " : " << e.t() << endl;
        iter++;

        //mise a jour des courbes
        vpPoseVector ctw(cTw);
        plot.plot(0, 0, iter, e.sumSquare());
        plot.plot(1, iter, e);
        plot.plot(2, iter, v);
        plot.plot(3, iter, ctw);
        //mise a jour de l'image
        display(cam, I, x, xd);
    }

    // sauvegarde des courbes
    plot.saveData(0, "e.txt", "#");
    plot.saveData(1, "error.txt", "#");
    plot.saveData(2, "v.txt", "#");
    plot.saveData(3, "p.txt", "#");

    // sauvegarde de l'image finale
    {
        vpImage<vpRGBa> Irgb;
        vpDisplay::getImage(I, Irgb);
        vpImageIo::write(Irgb, "4pt.jpg");
    }
    cout << "Clicker sur l'image pour terminer" << endl;
    vpDisplay::getClick(I);
}


void computeError3D(const vpHomogeneousMatrix& cdTc, vpColVector& e)
{
    e.resize(6);

    vpTranslationVector translation;
    cdTc.extract(translation);

    vpThetaUVector thetau;
    cdTc.extract(thetau);

    e[0]=translation[0];
    e[1]=translation[1];
    e[2]=translation[2];
    e[3]=thetau[0];
    e[4]=thetau[1];
    e[5]=thetau[2];
}

void preProduit(const vpColVector&u , vpMatrix& kU){
    kU.resize(3,3);
    kU=0;
    kU[0][1]=-u[2];
    kU[0][2]=u[1];

    kU[1][0]=u[2];
    kU[1][2]=-u[0];

    kU[2][0]=-u[1];
    kU[2][1]=u[0];
}

double sinc(double theta){
    if(theta==0)
        return 1;
    return sin(theta)/theta;
}


void computeInteractionMatrix3D(const vpHomogeneousMatrix& cdTc,vpMatrix &Lx)
{
    Lx=0;
    vpRotationMatrix rotation;
    cdTc.extract(rotation);

    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            Lx[i][j] = rotation[i][j];
        }
    }


    vpMatrix Lomega(3,3);
    Lomega.setIdentity();

    vpThetaUVector thetau;
    cdTc.extract(thetau);

    double theta;
    vpColVector u;
    thetau.extract(theta,u);

    vpMatrix Tu;
    preProduit(u,Tu);

    Lomega=Lomega+(theta/2.)*Tu+(1- sinc(theta)/(sinc(theta/2)*sinc(theta/2)))*Tu*Tu;

    std::cout << Lomega << std::endl;

    for (size_t i = 3; i < 6; i++)
    {
        for (size_t j = 3; j < 6; j++)
        {
            
            Lx[i][j] = Lomega[i-3][j-3];
        }
    }
}


void tp3DVisualServoing()
{

    vpTRACE("begin");

    vpPlot plot(4, 700, 700, 100, 200, "Curves...");

    char title[40];
    strncpy(title, "||e||", 40);
    plot.setTitle(0, title);
    plot.initGraph(0, 1);

    strncpy(title, "x-xd", 40);
    plot.setTitle(1, title);
    plot.initGraph(1, 6);

    strncpy(title, "camera velocity", 40);
    plot.setTitle(2, title);
    plot.initGraph(2, 6);

    strncpy(title, "Camera position", 40);
    plot.setTitle(3, title);
    plot.initGraph(3, 6);

    //Definition de la scene
    vpHomogeneousMatrix cTw(0, 0, 1.3, 0, 0, 0);
    //vpHomogeneousMatrix cTw(-0.2, -0.1, 1.3,vpMath::rad(10), vpMath::rad(20), vpMath::rad(30));
    //vpHomogeneousMatrix cTw(0, 0, 1, 0, 0, vpMath::rad(90));
    //vpHomogeneousMatrix cTw(0, 0, 1, 0, 0, vpMath::rad(180));

    vpHomogeneousMatrix cdTw(0, 0, 1, 0, 0, 0);

    vpColVector e(6);
    e=1;

    vpMatrix Lx(6, 6);

    vpColVector v(6);
    double lambda = 0.1;
    int iter = 0;
    while (fabs(e.sumSquare()) > 1e-6)
    {
        vpHomogeneousMatrix cdTc = cdTw*cTw.inverse();

        // Calcul de l'erreur
        computeError3D(cdTc,e) ;

        // Calcul de la matrice d'interaction
        computeInteractionMatrix3D(cdTc,Lx) ;

        //        Calcul de la loi de commande
        v = -lambda * Lx.pseudoInverse() * e;

        // Mis à jour de la position de la camera
        cTw = vpExponentialMap::direct(v).inverse() * cTw;

        cout << "iter " << iter << " : " << e.t() << endl;

        iter++;

        //mis a jour de courbes
        vpPoseVector crw(cTw);
        plot.plot(0, 0, iter, e.sumSquare());
        plot.plot(1, iter, e);
        plot.plot(2, iter, v);
        plot.plot(3, iter, crw);
    }

    // sauvegarde des courbes
    plot.saveData(0, "e.txt", "#");
    plot.saveData(1, "error.txt", "#");
    plot.saveData(2, "v.txt", "#");
    plot.saveData(3, "p.txt", "#");

    int a;
    cin >> a;
}

void tp2DVisualServoingFourPointMvt()
{

    //-------------------------------------------------------------
    // Mise en oeuvre des courbes

    vpPlot plot(4, 700, 700, 100, 200, "Curves...");

    char title[40];
    strncpy(title, "||e||", 40);
    plot.setTitle(0, title);
    plot.initGraph(0, 1);

    strncpy(title, "x-xd", 40);
    plot.setTitle(1, title);
    plot.initGraph(1, 8);

    strncpy(title, "camera velocity", 40);
    plot.setTitle(2, title);
    plot.initGraph(2, 6);

    strncpy(title, "camera position", 40);
    plot.setTitle(3, title);
    plot.initGraph(3, 6);

    //-------------------------------------------------------------
    // Affichage des images
    vpImage<unsigned char> I(400, 600);
    vpDisplayX d;
    d.init(I);
    vpDisplay::display(I);
    vpCameraParameters cam(400, 400, 300, 200);

    //-------------------------------------------------------------

    // Origine de la caméra

    vpHomogeneousMatrix cTw(0, 0, 1.3, 0, 0, 0);
    //vpHomogeneousMatrix cTw(-0.2, -0.1, 1.3,vpMath::rad(10), vpMath::rad(20), vpMath::rad(30));
    //vpHomogeneousMatrix cTw(0.2, 0.1, 1.3, 0, 0, vpMath::rad(5));
    //vpHomogeneousMatrix cTw(0, 0, 1, 0, 0, vpMath::rad(45));
    //vpHomogeneousMatrix cTw(0, 0, 1, 0, 0, vpMath::rad(90));
    //vpHomogeneousMatrix cTw(0, 0, 1, 0, 0, vpMath::rad(180));

    //-------------------------------------------------------------

    int nb_point=4;
    // Point cible dans le monde 3D

    double M = 0.40;

    std::srand(std::time(nullptr));

    vpColVector wX[nb_point];
    for (int i = 0; i < nb_point; i++){
        wX[i].resize(3);
        wX[i][0] = ((double)rand() / RAND_MAX)*2*M-M;
        wX[i][1] = ((double)rand() / RAND_MAX)*2*M-M;
        wX[i][2] = ((double)rand() / RAND_MAX)*2*M-M;
    }
    

    wX[0][0] = -M;
    wX[0][1] = -M;
    wX[0][2] = 0;
    wX[1][0] = M;
    wX[1][1] = -M;
    wX[1][2] = 0;
    wX[2][0] = M;
    wX[2][1] = M;
    wX[2][2] = 0;
    wX[3][0] = -M;
    wX[3][1] = M;
    wX[3][2] = 0;


    //-------------------------------------------------------------

    // Point à viser dans l'écran

    // position finale
    vpColVector xd(2*nb_point);
    vpHomogeneousMatrix cdTw(0, 0, 1, 0, 0, 0);

    //initialisation de la position désire des points dans l'image en fonction de cdTw
    for (size_t i = 0; i < nb_point; i++)
    {
        vpColVector cdX;
        changeFrame(wX[i], cdTw, cdX);
        vpColVector _xd;
        project(cdX, _xd);
        xd[i * 2 + 0] = _xd[0];
        xd[i * 2 + 1] = _xd[1];
    }
    
    //-------------------------------------------------------------

    vpColVector e(2*nb_point);
    e = 1;

    vpColVector x(2*nb_point);

    vpColVector v(6);
    double lambda = 0.1;
    int iter = 0;

    vpColVector v_obj(3);
    v_obj[0]=0;
    v_obj[1]=0;
    v_obj[2]=0.4;

    while (fabs(e.sumSquare()) > 1e-16)
    {
        std::cout << "----------------------------" << std::endl;
        std::cout << "Itération : " << iter << std::endl;

        for (int i = 0; i < nb_point; i++)
            wX[i]+=v_obj;

        std::cout << "  Déplacement de la cible"<< std::endl;
        vpMatrix Lx(2*nb_point, 6);

        // calcul de la position des points dans l'image en fonction de cTw
        for (size_t i = 0; i < nb_point; i++)
        {
            vpColVector _cX;
            changeFrame(wX[i], cTw, _cX);

            vpColVector _x;
            project(_cX, _x);

            x[i * 2 + 0] = _x[0];
            x[i * 2 + 1] = _x[1];

            vpMatrix _Lx;
            computeInteractionMatrix(_cX, _x[0], _x[1], _Lx);

            for (size_t j = 0; j < 6; j++)
            {
                Lx[i * 2 + 0][j] = _Lx[0][j];
                Lx[i * 2 + 1][j] = _Lx[1][j];
            }
        }

        //calcul de l'erreur
        e = x - xd;
        //calcul de la loi de commande
        v = -lambda * Lx.pseudoInverse() * e;
        //v[0]+=v_obj[0];
        //v[1]+=v_obj[1];
        //v[2]+=v_obj[2];

        //mise a jour de la position de la camera
        cTw = vpExponentialMap::direct(v).inverse() * cTw;

        //cout << "  iter " << iter << " : " << e.t() << endl;
        iter++;

        //mise a jour des courbes
        vpPoseVector ctw(cTw);
        plot.plot(0, 0, iter, e.sumSquare());
        plot.plot(1, iter, e);
        plot.plot(2, iter, v);
        plot.plot(3, iter, ctw);
        //mise a jour de l'image
        display(cam, I, x, xd);
    }

    // sauvegarde des courbes
    plot.saveData(0, "e.txt", "#");
    plot.saveData(1, "error.txt", "#");
    plot.saveData(2, "v.txt", "#");
    plot.saveData(3, "p.txt", "#");

    // sauvegarde de l'image finale
    {
        vpImage<vpRGBa> Irgb;
        vpDisplay::getImage(I, Irgb);
        vpImageIo::write(Irgb, "4pt.jpg");
    }
    cout << "Clicker sur l'image pour terminer" << endl;
    vpDisplay::getClick(I);
}

int main(int argc, char **argv)
{
    std::cout << "----------------------------" << std::endl;
    std::cout << "\n    TP Asservissement Visuel\n" << std::endl;
    std::cout << "----------------------------" << std::endl;
    std::cout << "1. Asservissement 2D" << std::endl;
    std::cout << "2. Asservissement 2D à 4 points" << std::endl;
    std::cout << "3. Asservissement 3D" << std::endl;
    std::cout << "4. Poursuite de points" << std::endl;

    std::cout << "Entrer un numéro :" << std::endl;
    int choice;
    std::cin >> choice;
    switch (choice)
    {
        case 1:
            tp2DVisualServoingOnePoint() ;
            break;
        case 2:
            tp2DVisualServoingFourPoint();
            break;
        case 3:
            tp3DVisualServoing() ;
            break;
        case 4:
            tp2DVisualServoingFourPointMvt();
            break;
        default:
            break;
    }
}
