#include "mesh.h"
#include "ui_vals.h"
#include <stdio.h>

using namespace boost::numeric::ublas;

// http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?LU_Matrix_Inversion
/* Matrix inversion routine.
   Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
bool InvertMatrix (matrix<double>& input, matrix<double>& inverse) {
   typedef permutation_matrix<std::size_t> pmatrix;
   // create a working copy of the input
   matrix<double> A(input);
   // create a permutation matrix for the LU-factorization
   pmatrix pm(A.size1());
   // perform LU-factorization
   int res = lu_factorize(A,pm);
       if( res != 0 ) return false;
   // create identity matrix of "inverse"
   inverse.assign(identity_matrix<double>(A.size1()));
   // backsubstitute to get the inverse
   lu_substitute(A, pm, inverse);
   return true;
}

mesh::mesh(QObject *parent) :
    QObject(parent)
{

}

 void mesh::init(){
     dx = user->x_len / (user->i_bar - 1);
     dy = user->y_len / (user->j_bar - 1);

     // create the matrices
     P = matrix<double>(user->i_bar + 2, user->j_bar + 2);
     psi = matrix<double>(user->i_bar + 2, user->j_bar + 2);
     u = matrix<double>(user->i_bar + 2, user->j_bar + 2);
     v = matrix<double>(user->i_bar + 2, user->j_bar + 2);
     u_b = matrix<double>(user->i_bar + 2, user->j_bar + 2);
     v_b = matrix<double>(user->i_bar + 2, user->j_bar + 2);
     D = matrix<double>(user->i_bar + 2, user->j_bar + 2);

     streamline_vals = matrix<double>(user->streamlines,1);

     // set them all to have the initial values
     for(int i = 0; i < user->i_bar + 2; i++)
     for(int j = 0; j < user->j_bar + 2; j++)
     {
         P(i,j) = user->P0;
         u(i,j) = user->u0;
         v(i,j) = user->v0;
     }

     // set boundary conditions

    // unless something weird is going on vis a vis blocked inputs or outputs, ul = ur
     for(int j = 0; j < user->j_bar+2; j++)
     {
         u(0,j) = user->ul;
         //u(1,j) = user->ul;
         //u(user->i_bar, j) = user->ur;
         //u(user->i_bar + 1, j) = user->ur;
         //u(user->i_bar - 1, j) = user->ur;
     }

     //something weird *is* up! comment this out when using other boundaries
     /*for(int j = user->j_bar/3; j < user->j_bar - user->j_bar/3; j++)
     {
         u(0,j) = 0;
     }*/

     // should be all 0s or other assumptions will break
     for(int i = 0; i < user->i_bar+1; i++)
     {
         v(i,0) = user->vb;
         v(i, user->j_bar) = user->vt;
     }

     beta = 1.0 / (2.0*user->dt * (1.0/(dx*dx) + 1.0/(dy*dy)));
     D_test = 0.002 * user->ul / dx;


     //circle_px(user->mast_r, user->mast_x, user->mast_y, &mast_pixels);
     circle_fill(user->mast_r, user->mast_x, user->mast_y, &mast_fill);
     mesh_mast();

     // add perturbations (todo)
 }

 void mesh::step()
 {
     explicit_step();
     for(int i = 0; i < 100 && implicit_step() > D_test; i++);
 }

 void mesh::explicit_step()
 {
    // set the tangential boundary conditions

    for(int i = 0; i < user->i_bar + 2; i++)
    {
        u(i,0) = u(i,1);
        u(i,user->j_bar) = u(i,user->j_bar-1);
        u(i,user->j_bar+1) = u(i,user->j_bar);

        v(i, 0) = 0;
        v(i,user->j_bar) = 0;
        v(i,user->j_bar+1) = 0;

    }
    for(int j = 0; j < user->j_bar + 2; j++)
    {
        v(0,j) = v(1,j);
        v(user->i_bar, j) = v(user->i_bar-1 , j);
    }

    /*printf("*****\n");
    // doesnt set values for the fictitious zones
    for(int i = 0; i < user->i_bar+2; i++)
    {
        for(int j = 0; j < user->j_bar+2; j++)
        {
            printf("(%2.5lf,%2.5lf)", u(i,j), v(i,j));
        }
        printf("\n");
    }
    printf("\n");*/

     // doesnt set values for the fictitious zones
     for(int i = 1; i < user->i_bar; i++)
     {
     for(int j = 1; j < user->j_bar; j++)
     {
         u_b(i,j) = u(i,j) - user->dt * (
                     (u(i+1,j)*u(i+1,j) - u(i-1,j)*u(i-1,j))/(2*dx) // du2/dx
                   + (u(i,j+1)*v(i,j+1) - u(i,j-1)*v(i,j-1))/(2*dy) // d(uv)/dy
                   - user->anu * (u(i+1,j) - 2*u(i,j) + u(i-1,j))/(dx*dx) // d2u/dx2
                   - user->anu * (u(i,j+1) - 2*u(i,j) + u(i,j-1))/(dy*dy) //d2u/dy2
                     );
         v_b(i,j) = v(i,j) - user->dt * (
                     (v(i,j+1)*v(i,j+1) - v(i,j-1)*v(i,j-1))/(2*dy) // dv2/dy
                   + (v(i+1,j)*u(i+1,j) - v(i-1,j)*u(i-1,j))/(2*dx) // d(uv)/dx
                   - user->anu * (v(i,j+1) - 2*v(i,j) + v(i,j-1))/(dy*dy) // d2v/dy2
                   - user->anu * (v(i+1,j) - 2*v(i,j) + v(i-1,j))/(dx*dx) //d2v/dx2
                     );
         //printf("(%2.5lf,%2.5lf)", u_b(i,j), v_b(i,j));
     }
     //printf("\n");
     }
     // printf("\n");
 }

 double mesh::implicit_step()
 {
    D_max = -1;

    // calculate D
    for(int i = 1; i < user->i_bar + 1; i++)
    for(int j = 1; j < user->j_bar + 1; j++)
    {
       D(i,j) = (u(i+1,j) - u(i-1,j))/(2*dx) + (v(i,j+1) - v(i,j-1))/(2*dy); // switch this to a centered difference?
       if(abs(D(i,j)) > D_max)
       {
           D_max = abs(D(i,j));
       }
    }
    for(int i = 0; i < user->i_bar+2; i++)
    {
       D(i,0) = D(i,1);
       D(i,user->j_bar) = D(i, user->j_bar - 1);
       D(i,user->j_bar+1) = D(i, user->j_bar);
    }
    for(int j = 0; j < user->j_bar+2; j++)
    {
       D(0,j) = D(1,j);
       D(user->i_bar,j) = D(user->i_bar - 1, j);
       D(user->i_bar+1,j) = D(user->i_bar, j);
    }
   // printf("D:\n");
    for(int i = 0; i < user->i_bar+2; i++)
    {
        for(int j = 0; j < user->j_bar+2; j++)
        {
            //printf("%2.5lf ", D(i,j));
        }
        //printf("\n");
    }
    //printf("\n");

    // calculate P
    for(int i = 1; i < user->i_bar+1; i++)
    for(int j = 1; j < user->j_bar+1; j++)
    {
       P(i,j) = P(i,j) - beta * D(i,j);
    }

    for(int i = 0; i < user->i_bar+1; i++)
    {
       P(i,0) = P(i,1);
       P(i,user->j_bar) = P(i, user->j_bar-1);
       P(i,user->j_bar+1) = P(i, user->j_bar);
    }
    for(int j = 0; j < user->j_bar+1; j++)
    {
       P(0,j) = P(1,j);
       P(user->i_bar,j) = P(user->i_bar - 1, j);
       P(user->i_bar+1,j) = P(user->i_bar, j);
    }

    //printf("P:\n");

    for(int i = 0; i < user->i_bar+1; i++)
    {
        for(int j = 0; j < user->j_bar+1; j++)
        {
           // printf("%2.5lf ", P(i,j));
        }
        //printf("\n");
    }
   // printf("\n");

    // calculate u
    for(int i = 1; i < user->i_bar+1; i++)
    for(int j = 1; j < user->j_bar+1; j++)
    {
        u(i,j) = u_b(i,j) - (user->dt/(2*dx)) * (P(i+1,j) - P(i - 1, j)); // switch to a centered difference?
    }

    // calculate v
    for(int i = 1; i < user->i_bar+1; i++)
    for(int j = 1; j < user->j_bar+1; j++)
    {
        v(i,j) = v_b(i,j) - (user->dt/(2*dy)) * (P(i,j+1) - P(i, j - 1)); // switch to a centered difference?
    }

    // mesh obstacles

    mesh_mast();




    // set boundary conditions
    /*sum_exit = 0;
    for(int j = 0; j < user->j_bar + 1; j++)
    {
        sum_exit += u(user->i_bar, j);
    }

    sum_entry = 0;
    for(int j = 1; j < user->j_bar; j++)
    {
        sum_entry += u(0, j);
    }*/
    //printf("\n");
    for(int j = 0; j < user->j_bar + 1; j++)
    {
        /*if(sum_exit != 0)
        {
            u(user->i_bar, j) = u(user->i_bar-1, j) * sum_entry / sum_exit;
            printf("%2.3lf\n", u(user->i_bar, j));
        }
        else // all zeroes, probably
        {*/
            u(user->i_bar, j) = u(user->i_bar - 1, j);
            //printf("%2.3lf\n", u(user->i_bar, j));
       // }
    }

    return D_max;
 }

 void mesh::mesh_mast()
 {
    for (std::vector<pixel_str>::iterator it = mast_pixels.begin() ; it != mast_pixels.end(); ++it)
    {
        mesh_pixel_bounds(*it);
    }
    // now set the pixel vels to 0;
    for (std::vector<pixel_str>::iterator it = mast_pixels.begin() ; it != mast_pixels.end(); ++it)
    {
        if((*it).ip >= 0 && (*it).jp >= 0 && (*it).ip <= user->i_bar+2 && (*it).jp <= user->j_bar + 2)
        {
            u((*it).ip, (*it).jp) = 0;
            v((*it).ip, (*it).jp) = 0;
        }
    }

    for (std::vector<pixel_str>::iterator it = mast_fill.begin() ; it != mast_fill.end(); ++it)
    {
        if((*it).ip >= 0 && (*it).jp >= 0 && (*it).ip <= user->i_bar+2 && (*it).jp <= user->j_bar + 2)
        {
            u((*it).ip, (*it).jp) = 0;
            v((*it).ip, (*it).jp) = 0;
        }
    }
 }

 void mesh::mesh_pixel_bounds(pixel_str px)
 {
    if(px.ip + 2 > user->i_bar + 1)
    {
        if(px.ip + 1 <= user->i_bar + 1)
        {
            v(px.ip+1, px.jp) = 0;
        }
    }
    else
    {
        v(px.ip+1, px.jp) = v(px.ip + 2, px.jp);
    }

    if(px.ip - 2 < 0)
    {
        if(px.ip - 1 >= 0)
        {
            v(px.ip-1, px.jp) = 0;
        }
    }
    else
    {
        v(px.ip-1, px.jp) = v(px.ip - 2, px.jp);
    }

    if(px.jp + 2 > user->j_bar + 1)
    {
        if(px.jp + 1 <= user->j_bar + 1)
        {
            u(px.ip, px.jp + 1) = 0;
        }
    }
    else
    {
        u(px.ip, px.jp + 1) = u(px.ip, px.jp + 2);
    }

    if(px.jp - 2 < 0)
    {
        if(px.jp - 1 >= 0)
        {
            u(px.ip, px.jp - 1) = 0;
        }
    }
    else
    {
        u(px.ip, px.jp - 1) = u(px.ip, px.jp - 2);
    }

 }

 // is actually an ellipse algorithm because dx = dy is typicaly false
void mesh::circle_px(double r, double xm, double ym, std::vector<pixel_str>*pixels)
 {
    pixels->clear();

    int x0 = (xm - r)/dx;
    int y0 = (ym - r)/dy;
    int x1 = (xm + r)/dx;
    int y1 = (ym + r)/dy;

    int a = abs(x1-x0), b = abs(y1-y0), b1 = b&1; /* values of diameter */
    long dx = 4*(1-a)*b*b, dy = 4*(b1+1)*a*a; /* error increment */
    long err = dx+dy+b1*a*a, e2; /* error of 1.step */

    if (x0 > x1) { x0 = x1; x1 += a; } /* if called with swapped points */
    if (y0 > y1) y0 = y1; /* .. exchange them */
    y0 += (b+1)/2; y1 = y0-b1;   /* starting pixel */
    a *= 8*a; b1 = 8*b*b;

    do {
        pixels->push_back({x1, y0}); /*   I. Quadrant */
        pixels->push_back({x0, y0}); /*  II. Quadrant */
        pixels->push_back({x0, y1}); /* III. Quadrant */
        pixels->push_back({x1, y1}); /*  IV. Quadrant */
        e2 = 2*err;
        if (e2 <= dy) { y0++; y1--; err += dy += a; }  /* y step */
        if (e2 >= dx || 2*err > dy) { x0++; x1--; err += dx += b1; } /* x step */
    } while (x0 <= x1);

    while (y0-y1 < b) {  /* too early stop of flat ellipses a=1 */
        pixels->push_back({x0-1, y0}); /* -> finish tip of ellipse */
        pixels->push_back({x1+1, y0++});
        pixels->push_back({x0-1, y1});
        pixels->push_back({x1+1, y1--});
    }
 }

void mesh::row(int x, int y, int width, std::vector<pixel_str>*pixels)
{
    for(int i = x; i < x + width; i++)
        pixels->push_back({i,y});
}

//http://enchantia.com/graphapp/doc/tech/ellipses.html
void mesh::circle_fill(double r, double xc, double yc, std::vector<pixel_str>*pixels)
{
     pixels->clear();

    /* e(x,y) = b^2*x^2 + a^2*y^2 - a^2*b^2 */
    int a = round(r/dx);
    int b = round(r/dy);
    xc = round(xc / dx);
    yc = round(yc / dy);


    #define incx() x++, dxt += d2xt, t += dxt
    #define incy() y--, dyt += d2yt, t += dyt
    int x = 0, y = b;
    unsigned int width = 1;
    long a2 = (long)a*a, b2 = (long)b*b;
    long crit1 = -(a2/4 + a%2 + b2);
    long crit2 = -(b2/4 + b%2 + a2);
    long crit3 = -(b2/4 + b%2);
    long t = -a2*y; /* e(x+1/2,y-1/2) - (a^2+b^2)/4 */
    long dxt = 2*b2*x, dyt = -2*a2*y;
    long d2xt = 2*b2, d2yt = 2*a2;

    while (y>=0 && x<=a) {
        if (t + b2*x <= crit1 ||     /* e(x+1,y-1/2) <= 0 */
            t + a2*y <= crit3) {     /* e(x+1/2,y) <= 0 */
            incx();
            width += 2;
        }
        else if (t - a2*y > crit2) { /* e(x+1/2,y-1) > 0 */
            row(xc-x, yc-y, width, pixels);
            if (y!=0)
                row(xc-x, yc+y, width, pixels);
            incy();
        }
        else {
            row(xc-x, yc-y, width, pixels);
            if (y!=0)
                row(xc-x, yc+y, width, pixels);
            incx();
            incy();
            width += 2;
        }
    }
    if (b == 0)
        row(xc-a, yc, 2*a+1, pixels);
}




 void mesh::draw_streamlines_object()
 {

     psi_max = 0;

     for(int i = 0 ; i < user->i_bar + 1; i++)
     {
         psi(i,0) = 0;
     }

     for(int i = 0; i < user->i_bar + 1; i++)
     for(int j = 1; j < user->j_bar + 1; j++)
     {
         psi(i,j) = psi(i, j-1) + u(i,j) * dy;
         if(psi(i,j) > psi_max)
         {
             psi_max = psi(i,j);
         }
     }
 }

 void mesh::draw_streamlines_fluid()
 {
     psi_max = 0;

     for(int i = 0 ; i < user->i_bar + 1; i++)
     {
         psi(i,0) = 0;
     }

     for(int i = 0; i < user->i_bar + 1; i++)
     for(int j = 1; j < user->j_bar + 1; j++)
     {
         psi(i,j) = psi(i, j-1) + (u(i,j) - user->ur) * dy;
         if(psi(i,j) > psi_max)
         {
             psi_max = psi(i,j);

         }
     }
 }

 // leave this for later
 void mesh::draw_pressure_field()
 {
     double r = 0;
     double g = 0;
     double b = 0;

     glEnable(GL_POINT_SMOOTH);
     glPointSize(3);

     for(int i = 1; i < user->i_bar+1; i++)
     for(int j = 1; j < user->j_bar+1; j++)
     {
         //printf("%2.3lf\n",u(i,j));
         if(u(i,j) > user->P_max/2)
         {
            r = 2*(255.0 - 190)*(u(i,j) - user->P_max/2)/(user->P_max) + 190;
            if(r > 255) r = 255;
         }
         else if(u(i,j) > 0)
         {
             r = 2*(190.0 - 100)*(u(i,j))/(user->P_max) + 100;
             g = 2*(0 - 100)*(u(i,j))/(user->P_max) + 100;
             b = 2*(0 - 100)*(u(i,j))/(user->P_max) + 100;
         }
         else if(u(i,j) > user->P_min/2)
         {
             b = 2*(190.0 - 100)*(u(i,j))/(user->P_min) + 100;
             g = 2*(0 - 100)*(u(i,j))/(user->P_min) + 100;
             r = 2*(0 - 100)*(u(i,j))/(user->P_min) + 100;
         }
         else
         {
             b = 2*(255.0 - 190)*(u(i,j) - user->P_min/2)/(user->P_min) + 190;
             if(b > 255) b = 255;
         }
         glColor3f(r/255.0, g/255.0, b/255.0);

         glBegin(GL_POINTS);
         glVertex2f(i*dx, j*dy);
         glEnd();
     }
     glColor3f(0.0, 0.0, 0.0);
 }


// http://paulbourke.net/papers/conrec/conrec.c, simplified to get rid of some values and inputs
 /*
    Derivation from the fortran version of CONREC by Paul Bourke
    d               ! matrix of data to contour
    ilb,iub,jlb,jub ! index bounds of data matrix
    x               ! data matrix column coordinates
    y               ! data matrix row coordinates
    nc              ! number of contour levels
    z               ! contour levels in increasing order
 */
 void mesh::draw_contour()
 {
 #define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
 #define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])

    int ilb = 0;
    int iub = user->i_bar + 1;
    int jlb = 0;
    int jub = user->j_bar + 1;

    int m1,m2,m3,case_value;
    double dmin,dmax,x1=0,x2=0,y1=0,y2=0;
    int i,j,k,m;
    double h[5];
    int sh[5];
    double xh[5],yh[5];
    int im[4] = {0,1,1,0},jm[4]={0,0,1,1};
    int castab[3][3][3] = {
      { {0,0,8},{0,2,5},{7,6,9} },
      { {0,3,4},{1,3,1},{4,3,0} },
      { {9,6,7},{5,2,0},{8,0,0} }
    };
    double temp1,temp2;

    for(k = 0; k < user->streamlines; k++)
    {
        streamline_vals(k,0) = ((double)k) / user->streamlines * psi_max;
    }

    for (j=(jub-1);j>=jlb;j--) {
       for (i=ilb;i<=iub-1;i++) {
          temp1 = fmin(psi(i,j),psi(i,j+1));
          temp2 = fmin(psi(i+1,j),psi(i+1,j+1));
          dmin  = fmin(temp1,temp2);
          temp1 = fmax(psi(i,j),psi(i,j+1));
          temp2 = fmax(psi(i+1,j),psi(i+1,j+1));
          dmax  = fmax(temp1,temp2);
          if (dmax < streamline_vals(0,0) || dmin > streamline_vals(user->streamlines - 1,0))
             continue;
          for (k=0;k<user->streamlines; k++) {
             if (streamline_vals(k,0) < dmin || streamline_vals(k,0) > dmax)
                continue;
             for (m=4;m>=0;m--) {
                if (m > 0) {
                   h[m]  = psi(i+im[m-1],j+jm[m-1])-streamline_vals(k,0);
                   xh[m] = i+im[m-1];
                   yh[m] = j+jm[m-1];
                } else {
                   h[0]  = 0.25 * (h[1]+h[2]+h[3]+h[4]);
                   xh[0] = 0.50 * (i+i+1);
                   yh[0] = 0.50 * (j+j+1);
                }
                if (h[m] > 0.0)
                   sh[m] = 1;
                else if (h[m] < 0.0)
                   sh[m] = -1;
                else
                   sh[m] = 0;
             }

             /*
                Note: at this stage the relative heights of the corners and the
                centre are in the h array, and the corresponding coordinates are
                in the xh and yh arrays. The centre of the box is indexed by 0
                and the 4 corners by 1 to 4 as shown below.
                Each triangle is then indexed by the parameter m, and the 3
                vertices of each triangle are indexed by parameters m1,m2,and m3.
                It is assumed that the centre of the box is always vertex 2
                though this isimportant only when all 3 vertices lie exactly on
                the same contour level, in which case only the side of the box
                is drawn.
                   vertex 4 +-------------------+ vertex 3
                            | \               / |
                            |   \    m-3    /   |
                            |     \       /     |
                            |       \   /       |
                            |  m=2    X   m=2   |       the centre is vertex 0
                            |       /   \       |
                            |     /       \     |
                            |   /    m=1    \   |
                            | /               \ |
                   vertex 1 +-------------------+ vertex 2
             */
             /* Scan each triangle in the box */
             for (m=1;m<=4;m++) {
                m1 = m;
                m2 = 0;
                if (m != 4)
                   m3 = m + 1;
                else
                   m3 = 1;
                if ((case_value = castab[sh[m1]+1][sh[m2]+1][sh[m3]+1]) == 0)
                   continue;
                switch (case_value) {
                case 1: /* Line between vertices 1 and 2 */
                    x1 = xh[m1];
                    y1 = yh[m1];
                    x2 = xh[m2];
                    y2 = yh[m2];
                    break;
                case 2: /* Line between vertices 2 and 3 */
                    x1 = xh[m2];
                    y1 = yh[m2];
                    x2 = xh[m3];
                    y2 = yh[m3];
                    break;
                case 3: /* Line between vertices 3 and 1 */
                    x1 = xh[m3];
                    y1 = yh[m3];
                    x2 = xh[m1];
                    y2 = yh[m1];
                    break;
                case 4: /* Line between vertex 1 and side 2-3 */
                    x1 = xh[m1];
                    y1 = yh[m1];
                    x2 = xsect(m2,m3);
                    y2 = ysect(m2,m3);
                    break;
                case 5: /* Line between vertex 2 and side 3-1 */
                    x1 = xh[m2];
                    y1 = yh[m2];
                    x2 = xsect(m3,m1);
                    y2 = ysect(m3,m1);
                    break;
                case 6: /* Line between vertex 3 and side 1-2 */
                    x1 = xh[m3];
                    y1 = yh[m3];
                    x2 = xsect(m1,m2);
                    y2 = ysect(m1,m2);
                    break;
                case 7: /* Line between sides 1-2 and 2-3 */
                    x1 = xsect(m1,m2);
                    y1 = ysect(m1,m2);
                    x2 = xsect(m2,m3);
                    y2 = ysect(m2,m3);
                    break;
                case 8: /* Line between sides 2-3 and 3-1 */
                    x1 = xsect(m2,m3);
                    y1 = ysect(m2,m3);
                    x2 = xsect(m3,m1);
                    y2 = ysect(m3,m1);
                    break;
                case 9: /* Line between sides 3-1 and 1-2 */
                    x1 = xsect(m3,m1);
                    y1 = ysect(m3,m1);
                    x2 = xsect(m1,m2);
                    y2 = ysect(m1,m2);
                    break;
                default:
                    break;
                }

                /* Finally draw the line */
                // we don't really need anything fancy for this, just to plot the new line segment.
                // Alltogether, it should workout without needing to make a polyline to follow
                glBegin(GL_LINES);
                    glVertex2f(x1*dx, y1*dy);
                    glVertex2f(x2*dx, y2*dy);
                glEnd();

             } /* m */
          } /* k - contour */
       } /* i */
    } /* j */
 }

 void mesh::draw_mast()
 {
     float cx = user->mast_x;
     float cy = user->mast_y;
     float r = user->mast_r;
     int num_segments = 360;
     float theta = 2 * 3.1415926 / float(num_segments);
     float tangetial_factor = tanf(theta);//calculate the tangential factor

     float radial_factor = cosf(theta);//calculate the radial factor

     float x = r;//we start at angle = 0

     float y = 0;

     glBegin(GL_POLYGON);
     for(int ii = 0; ii < num_segments; ii++)
     {
         glVertex2f(x + cx, y + cy);//output vertex

         //calculate the tangential vector
         //remember, the radial vector is (x, y)
         //to get the tangential vector we flip those coordinates and negate one of them

         float tx = -y;
         float ty = x;

         //add the tangential vector

         x += tx * tangetial_factor;
         y += ty * tangetial_factor;

         //correct using the radial factor

         x *= radial_factor;
         y *= radial_factor;
     }
     glEnd();
    }
