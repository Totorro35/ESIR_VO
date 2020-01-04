#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <visp3/core/vpImage.h>
#include <visp3/io/vpImageIo.h>
#include "EigenFacesDB.h"
#include <omp.h>

void  writeEigenFaceInImage(int m_w, int m_h, const std::string& path, vpColVector &v)
{
    double maxV = v.getMaxValue();
    double minV = v.getMinValue();
    vpImage<unsigned char> I(m_h,m_w);
    for(int i=0; i<m_h; i++)
        for(int j = 0; j<m_w; j++)
        {
            I[i][j] = ((v[i*m_w+j]-minV)/(maxV-minV))*255;
        }
    vpImageIo::write(I, path);
}

void buildMeanImage(const vpMatrix& dataBase,vpColVector& m_vMean){
    int size = dataBase.getRows();
    m_vMean = vpColVector(size);
    for (int j = 0; j < dataBase.getCols(); j++)
    {
        m_vMean+=(dataBase.getCol(j)/dataBase.getCols());
    }
}

void buildBDFaces(const std::vector<std::string>& paths, int maxEigenFace)
{
    int m_maxEigenFace = std::min(maxEigenFace, (int)paths.size());

    int K = m_maxEigenFace/2;
    int nbTirage = 4;
    
    // On calcul les attibuts de l'image
    vpImage<unsigned char> I;
    vpImageIo::read(I,*paths.begin());
    int m_w = I.getWidth();
    int  m_h = I.getHeight();
    int m_size = m_h*m_w;
    
    std::cout << " * Caracteristique de l'images : " << m_h << "x" << m_w << std::endl;
    std::cout << " * Nombre d'image de la base : " << paths.size()<< std::endl;
    std::cout << " * Nombre de U : " << m_maxEigenFace << std::endl;


    // Creation du vpColVector pour le mean face
    std::cout << "[INFO] Loading Reference Base" << std::endl;
    vpMatrix dataBase(m_size,m_maxEigenFace);
    int cpt=0;
    for(auto path : paths)
    {
        if(cpt>m_maxEigenFace)
            break;

        vpImageIo::read(I,path);
        for (int i = 0; i < m_h; i++)
        {
            for (int j = 0; j < m_w; j++)
            {
                dataBase[i*m_w+j][cpt] = I[i][j]/255.;
            }
        }
        cpt++;
    }

    // Creation du vpColVector pour le mean face
    std::cout << "[INFO] Loading Test Base" << std::endl;
    vpMatrix dataBase_Test(m_size,40);
    for (size_t idx = 0; idx < 40; idx++)
    {
        std::string path = "../Donnees/att_faces/s"+std::to_string(idx+1)+"/10.pgm";
        vpImageIo::read(I,path);
        for (int i = 0; i < m_h; i++)
        {
            for (int j = 0; j < m_w; j++)
            {
                dataBase_Test[i*m_w+j][idx] = I[i][j]/255.;
            }
        }
    }
    
    // Creation du vpColVector pour le mean face
    std::cout << "[INFO] Creation Mean images ...." << std::endl;
    vpColVector m_vMean ;
    buildMeanImage(dataBase,m_vMean);

    writeEigenFaceInImage(m_w,m_h,"../result/mean.jpg",m_vMean);
    srand(time(0));  
    for (int i = 0; i < 4; i++)
    {
        int elem = rand()%m_maxEigenFace;
        vpColVector face = dataBase.getCol(elem);
        std::string save_path = "../result/"+std::to_string(i)+".jpg";
        writeEigenFaceInImage(m_w,m_h,save_path,face);

        save_path = "../result/"+std::to_string(i)+"_centre.jpg";
        vpColVector centre = face-m_vMean;
        writeEigenFaceInImage(m_w,m_h,save_path,centre);
    }
    
    // Calcul de la matrice A
    std::cout << "[INFO] Calcul de A ... " << std::endl;
    vpMatrix A = vpMatrix(m_size,m_maxEigenFace);
    for (int j = 0; j < m_maxEigenFace; j++)
    {
        for(int i=0;i<m_size;i++)
        {
            A[i][j]=dataBase[i][j]-m_vMean[i];
        }
    }

    // SVD
    std::cout << "[INFO] Calcul de U ... " << std::endl;
    vpMatrix v;
    vpColVector eigenValue;
    vpMatrix eigenVector=A;
    eigenVector.svd(eigenValue,v);

    std::cout << "[INFO] Sauvegarde des 4 premières eigenfaces " << std::endl;
    for (int j = 0; j < 4; j++)
    {
        vpColVector eigenface = eigenVector.getCol(j);
        std::string save_path = "../result/eigen_face_"+std::to_string(eigenValue[j])+".jpg";
        writeEigenFaceInImage(m_w,m_h,save_path,eigenface);
    }

    std::cout << "[INFO] Reprojection d'une image de reference pour K=50,150 et 360 " << std::endl;
    int K_values[3] = {50,150,360};

    int elem = rand()%m_maxEigenFace; // Choix d'un visage aléatoire
    vpColVector face = dataBase.getCol(elem);

    for (size_t i = 0; i < 3; i++)
    {
        K = K_values[i];
        vpColVector W(K);
        for (size_t k = 0; k < K; k++)
        {
            W[k]=eigenVector.getCol(k)*(face-m_vMean);
        }

        vpColVector reproj_I = vpColVector(m_size);
        for (size_t k = 0; k < K; k++)
        {
            reproj_I+=W[k]*eigenVector.getCol(k);
        }
        reproj_I+=m_vMean;

        double error = (face-reproj_I).infinityNorm();
        std::cout << "K="<<K<<" ; Error=" << error << std::endl;

        std::string save_path = "../result/reproj_ref_"+std::to_string(K)+".jpg";
        writeEigenFaceInImage(m_w,m_h,save_path,reproj_I);
    }

    std::cout << "[INFO] Reprojection d'une image de test pour K=50,150 et 360 " << std::endl;
    elem = elem/9; // Visage de la même personne pour le test
    face = dataBase_Test.getCol(elem);

    for (size_t i = 0; i < 3; i++)
    {
        K = K_values[i];
        vpColVector W(K);
        for (size_t k = 0; k < K; k++)
        {
            W[k]=eigenVector.getCol(k)*(face-m_vMean);
        }

        vpColVector reproj_I = vpColVector(m_size);
        for (size_t k = 0; k < K; k++)
        {
            reproj_I+=W[k]*eigenVector.getCol(k);
        }
        reproj_I+=m_vMean;

        double error = (face-reproj_I).infinityNorm();
        std::cout << "K="<<K<<" ; Error=" << error << std::endl;

        std::string save_path = "../result/reproj_test_"+std::to_string(K)+".jpg";
        writeEigenFaceInImage(m_w,m_h,save_path,reproj_I);
    }


    std::cout << "[INFO] Calcul Normalized Cumulative EigenValue " << std::endl;
    std::cout << "Display result ? (Y/N)"<<std::endl;
    std::string display;
    std::cin >> display;

    if(display=="Y"){
        double sum_eigenvalue = eigenValue.sum();
        vpColVector cumule_eigenvalue = eigenValue;
        for (size_t i = 0; i < cumule_eigenvalue.size(); i++)
        {
            if(i==0){
                cumule_eigenvalue[i] = eigenValue[i] / sum_eigenvalue;
                continue;
            }else{
                cumule_eigenvalue[i] = cumule_eigenvalue[i-1] + eigenValue[i] / sum_eigenvalue;
            }
            std::cout << cumule_eigenvalue[i] << std::endl;
        }
    }


    std::cout << "[INFO] Error reconstruction " << std::endl;
    std::cout << "Display result ? (Y/N)"<<std::endl;
    std::cin >> display;

    if(display=="Y"){

    for (size_t nbEigen = 0; nbEigen < m_maxEigenFace; nbEigen++)
    {
        double error=0;

#pragma omp parallel for
        for (size_t elem = 0; elem < m_maxEigenFace; elem++)
        {
            vpColVector face = dataBase.getCol(elem);

            vpColVector W(nbEigen);
            for (size_t k = 0; k < nbEigen; k++)
            {
                W[k]=eigenVector.getCol(k)*(face-m_vMean);
            }

            vpColVector reproj_I = vpColVector(m_size);
            for (size_t k = 0; k < nbEigen; k++)
            {
                reproj_I+=W[k]*eigenVector.getCol(k);
            }
            reproj_I+=m_vMean;

#pragma omp critical
            error += (face-reproj_I).infinityNorm();
        }

        error/=m_maxEigenFace;
        std::cout << error << std::endl;
    }
    }

    std::cout << "[INFO] Compute Error Matrix for K=20 " << std::endl;
    std::cout << "Display result ? (Y/N)"<<std::endl;
    std::cin >> display;

    if(display=="Y"){
    int K_value=20;

    vpMatrix ErrorMatrix(m_maxEigenFace,m_maxEigenFace);

#pragma omp parallel for
    for (size_t i = 0; i < m_maxEigenFace; i++)
    {
        vpColVector face_i = dataBase.getCol(i);

        vpColVector W_i(K_value);
        for (size_t k = 0; k < K_value; k++)
        {
            W_i[k]=eigenVector.getCol(k)*(face_i-m_vMean);
        }

        for (size_t j = 0; j < m_maxEigenFace; j++)
        {
            vpColVector face_j = dataBase.getCol(j);

            vpColVector W_j(K_value);
            for (size_t k = 0; k < K_value; k++)
            {
                W_j[k]=eigenVector.getCol(k)*(face_j-m_vMean);
            }
        
            double error = (W_i - W_j).infinityNorm();

#pragma omp critical
            ErrorMatrix[i][j] = error;
        }
    }

    std::ofstream myfile;
    myfile.open ("../result/matrixError.txt");
    myfile<<ErrorMatrix;
    myfile.close();

    double min_same = 50000;
    double max_same = 0;
    double min_diff = 50000;
    double max_diff = 0;

    for (size_t i = 0; i < m_maxEigenFace; i++)
    {
        for (size_t j = 0; j < m_maxEigenFace; j++)
        {
            double value = ErrorMatrix[i][j];
            if(i/9==j/9){
                min_same=std::min(value,min_same);
                max_same=std::max(value,max_same);
            }else{
                min_diff=std::min(value,min_diff);
                max_diff=std::max(value,max_diff);
            }
        }
    }

    std::cout << "Min Same : "<<min_same << std::endl;
    std::cout << "Max Same : "<<max_same << std::endl;
    std::cout << "Min Diff : "<<min_diff << std::endl;
    std::cout << "Max Diff : "<<max_diff << std::endl;

    double theta = (max_same+min_diff)/2;
    std::cout << "Theta : "<<theta<<std::endl;
    }

    std::cout << "[INFO] Reconnaissance facial en fonction de K " << std::endl;
    std::cout << "Display result ? (Y/N)"<<std::endl;
    std::cin >> display;

    if(display=="Y"){
    double theta_seuil = 8.5;
    std::cout<<"Theta : "<<theta_seuil<<std::endl;

    for (size_t K_value = 0; K_value < m_maxEigenFace; K_value++)
    {
        std::cout << "  K : "<<K_value<< std::endl;
        int nbReconnu=0;

#pragma omp parallel for
        for (size_t i = 0; i < 40; i++)
        {
            vpColVector face_i = dataBase_Test.getCol(i);

            vpColVector W_i(K_value);
            for (size_t k = 0; k < K_value; k++)
            {
                W_i[k]=eigenVector.getCol(k)*(face_i-m_vMean);
            }

            for (size_t j = 0; j < m_maxEigenFace; j++)
            {
                vpColVector face_j = dataBase.getCol(j);

                vpColVector W_j(K_value);
                for (size_t k = 0; k < K_value; k++)
                {
                    W_j[k]=eigenVector.getCol(k)*(face_j-m_vMean);
                }
        
                double error = (W_i - W_j).infinityNorm();
                if(error<theta_seuil){
#pragma omp critical
                    nbReconnu++;
                    break;
                }
            }
        }
        std::cout << "    Nb Reconnu : "<<nbReconnu<< std::endl;
    }
    }

    std::cout << "[INFO] Fin du programme" << std::endl;
}

std::vector<std::string> buildPathImagesAttFaces()
{
	std::vector<std::string> v;
	for(int nbDir=1; nbDir<=40; nbDir++)
		for(int nbImage=1; nbImage<10;nbImage++)
		{
			std::ostringstream ss;
			ss << "../Donnees/att_faces/s" << nbDir << "/" << nbImage << ".pgm";
			v.push_back(ss.str());
		}
	return v;
}

int main()
{
	std::cout << "[INFO] Construction du path ..." << std::endl;
	std::vector<std::string> paths = buildPathImagesAttFaces();
	
    buildBDFaces(paths,400);
	
	return 0;
}
