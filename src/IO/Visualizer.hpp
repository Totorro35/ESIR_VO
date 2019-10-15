/**
 * @file Visualizer.hpp
 * @author Thomas LEMETAYER (thomas.lemetayer.35@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2019-09-27
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#pragma once

#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp3/gui/vpDisplayX.h>

/**
 * @brief Namespace IO
 * 
 */
namespace IO
{

/**
 * @brief Function to show an Image
 * 
 * @param title 
 * @param I the Image in RGB
 */
void show(const std::string &title, vpImage<vpRGBa> &I)
{
    vpDisplayX d(I, vpDisplay::SCALE_AUTO);
    vpDisplay::setTitle(I, title);
    vpDisplay::display(I);
    vpDisplay::flush(I);
    vpDisplay::getClick(I);
}

/**
 * @brief Main Function to test show()
 * 
 */
void MainVisu()
{
    vpImage<vpRGBa> I;
    vpImageIo::read(I, "../data/test.png");
    IO::show("C&V", I);
}

} // namespace IO