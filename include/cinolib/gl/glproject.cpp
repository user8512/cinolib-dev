/********************************************************************************
*  This file is part of CinoLib                                                 *
*  Copyright(C) 2022: Marco Livesu                                              *
*                                                                               *
*  The MIT License                                                              *
*                                                                               *
*  Permission is hereby granted, free of charge, to any person obtaining a      *
*  copy of this software and associated documentation files (the "Software"),   *
*  to deal in the Software without restriction, including without limitation    *
*  the rights to use, copy, modify, merge, publish, distribute, sublicense,     *
*  and/or sell copies of the Software, and to permit persons to whom the        *
*  Software is furnished to do so, subject to the following conditions:         *
*                                                                               *
*  The above copyright notice and this permission notice shall be included in   *
*  all copies or substantial portions of the Software.                          *
*                                                                               *
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR   *
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
*  FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE *
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      *
*  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS *
*  IN THE SOFTWARE.                                                             *
*                                                                               *
*  Author(s):                                                                   *
*                                                                               *
*     Marco Livesu (marco.livesu@gmail.com)                                     *
*     http://pers.ge.imati.cnr.it/livesu/                                       *
*                                                                               *
*     Italian National Research Council (CNR)                                   *
*     Institute for Applied Mathematics and Information Technologies (IMATI)    *
*     Via de Marini, 6                                                          *
*     16149 Genoa,                                                              *
*     Italy                                                                     *
*********************************************************************************/
#include <cinolib/gl/glproject.h>

namespace cinolib
{

// this is equivalent to gluProject
// https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluProject.xml
// https://www.khronos.org/opengl/wiki/GluProject_and_gluUnProject_code
// https://stackoverflow.com/questions/54385727/gluunproject-giving-incorrect-values-on-macretina-display-but-works-correctl
//
template<class T>
CINO_INLINE
void gl_project(const mat<3,1,T>   & p_in,       // 3d point in object space
                const mat<4,4,T>     modelview,  //
                const mat<4,4,T>     projection, //
                const mat<2,2,int>   viewport,   //
                      mat<2,1,T>   & p_out,      // 2d point in screen space
                      T            & depth)      // z-buffer
{
    mat<4,1,T> p   = p_in.add_coord(1); // add one (homogeneous) coordinate
    mat<4,1,T> res = (projection * modelview) * p;
    res /= res[3]; // convert back to cartesian coordinates

    p_out[0] = viewport[0] + (viewport[2]*(res[0]+1))/2.0;
    p_out[1] = viewport[1] + (viewport[3]*(res[1]+1))/2.0;
    depth    = (res[2]+1)/2.0;
}

}
