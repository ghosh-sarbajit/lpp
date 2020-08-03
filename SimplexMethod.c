#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define P 30
#define Q 30

typedef struct
{
    // p = n_rows, q = n_columns, mat[p x q] of augmented matrix
    int p,q;
    double matrix[P][Q];
} lpp;

static const double diff=1.0e-6;
int is_equal(double x,double y)
{
    return fabs(x-y) < diff;
}

//for printing a dashed line
void line_print(int num)
{
    for(int i=0;i<num;i++)
        putchar('-');
    printf("\n");
}

void print_problem_step(lpp *problem, const char *status)
{
    static int count=0;
    printf("\n%d. Problem %s:\n", ++count, status);
    line_print(75);
    // '-' is used to print in front
    printf("%-6s%5s", "col:", "b[i]");
    //to print varable name
    for(int j=1; j<problem->q; j++)
    {
        printf("    x%d ", j);
    }
    printf("\n" );
    //print data; p=n_rows, q=n_columns, matrix[pxq]
    for(int i=0;i<problem->p;i++)
    {
        if(i==0)
            printf("max:" );
        else
            printf("b%d: ", i);
        for(int j=0;j<problem->q; j++)
        {
            //if interger then interger ow double
            if (is_equal((int)problem->matrix[i][j], problem ->matrix[i][j]))
              printf(" %6d", (int)problem->matrix[i][j]);
            else
              printf(" %6.2lf", problem->matrix[i][j]);
        }
        printf("\n");
    }
    line_print(75);
}

void add_slack_vars(lpp *problem)
{
    for(int i=1; i<problem->p; i++)
    {
        for(int j=1; j<problem->p; j++)
        {
            problem->matrix[i][j + problem->q -1] = (i==j);
        }
    }
    problem->q += problem->p - 1;
}

// takes a col of identity mat, find the row containing 1.
// return -1 if col is not from an identity mat
int find_basis_vars(lpp *problem, int col)
{
    int xi=-1;
    for(int i=1; i< problem->p; i++)
    {
        if( is_equal(problem->matrix[i][col],1) )
        {
            if(xi == -1)
                xi = i;
            else
                return -1;
        }
        else if( !is_equal(problem->matrix[i][col],0) )
        {
            return -1;
        }
    }
    return xi;
}


void print_optimal_vec(lpp *problem, const char *status)
{
    int xi;
    printf("%s at\n", status);
    for(int j=1; j<problem->q;j++)
    {
        // will go through each column
        xi = find_basis_vars(problem, j);
        if(xi != -1)
            printf("x%d=%3.2lf ", j, problem->matrix[xi][0]);
        else
            printf("x%d=0, ", j);
    }
    printf("\n" );
}

void is_b_positive(lpp *problem)
{
    for(int i=1; i<problem->q; i++)
        assert(problem->matrix[i][0] >= 0);
}


// row operation
void pivot(lpp *problem, int row, int col)
{
    //row and col set by pivot
    double pivot;
    pivot = problem-> matrix[row][col];
    for(int j=0;j<problem->q;j++)
    {
        problem->matrix[row][j] /= pivot;
    }
    assert(is_equal(problem->matrix[row][col],1));

    for(int i=0; i<problem->p; i++)
    {
        //for other rows
        double k = problem-> matrix[i][col];
        if(i==row)
            continue;
        for(int j=0; j<problem->q;j++)
        {
            // r[i] = r[i] - k * r[row]
            problem->matrix[i][j] -= k * problem -> matrix[row][j];
        }
    }
}

// find pivot col index
// most negative column in matrix[0][1:Q]
int find_pivot_c_index(lpp *problem)
{
    int pivot_c = 1;
    double least_val = problem -> matrix[0][pivot_c];
    for(int j=1;j<problem->q; j++)
    {
        if (problem -> matrix[0][j] < least_val)
        {
            least_val = problem->matrix[0][j];
            pivot_c = j;
        }
    }
    printf("Most -ve col in 0_th row is col %d = %g\n",
                pivot_c, least_val);
    if(least_val >= 0)
    {
        // row[0] all positive
        // optimality reached !!!!
        return -1;
    }
    return pivot_c;
}

// find pivot row index
// row with min {theta}
// theta = col[0]/col[pivot]
int find_pivot_r_index(lpp *problem, int pivot_c)
{
    int pivot_r = 0;
    double min_theta = -1;
    printf("theta x[row_i,0]/x[row_i,%d] = [", pivot_c);
    for(int i=1; i< problem->p; i++)
    {
        double theta = problem->matrix[i][0]/ problem->matrix[i][pivot_c];
        printf("%3.2lf, ", theta);
        if( (theta > 0 && theta < min_theta) || min_theta <0)
        {
            min_theta = theta;
            pivot_r = i;
        }
    }
    printf("].\n" );
    // unbounded case
    if(min_theta == -1)
        return -1;
    printf("pivot at A[%d,%d], min positive ratio = %g in row %d\n",
            pivot_r, pivot_c, min_theta, pivot_r);
    return pivot_r;
}

void solve_by_simplex(lpp *problem)
{
  int loop=0;
  add_slack_vars(problem);
  is_b_positive(problem);
  //print augmented matrix after adding slack variables
  print_problem_step(problem,"slack variables added");
  while (++loop)
  {
      int pivot_c,pivot_r;
      pivot_c = find_pivot_c_index(problem);
      if(pivot_c < 0)
      {
        //row 0 all positive, optimality !!!!
        printf("opt value is at A[0,0]=%3.2lf\n", problem->matrix[0][0]);
        print_optimal_vec(problem,"optimal solution");
        break;
      }
      printf("Entering var x%d, into basis, pivot_c =%d\n", pivot_c,pivot_c);
      pivot_r = find_pivot_r_index(problem, pivot_c);
      if(pivot_r < 0)
      {
          printf("unbounded problem\n");
          break;
      }
      printf("Exiting var x%d, pivot_r =%d\n", pivot_r,pivot_r);
      pivot(problem,pivot_r,pivot_c);
      print_problem_step(problem,"after row operation");
      print_optimal_vec(problem,"BFS");
      if(loop > 30)
      {
          printf("Too many(%d) iterations\n", loop);
          break;
      }
  }
}


lpp problem  = { 4, 5, {              // problem size 4*5
    {  0.0 , -1.5 , -3.0 ,-1.0 , -4.0,   },  // Max: z = -1.5*x + 3*y + z + 4*w,
    { 20.0 ,  1.0 ,  1.0 , 1.0 ,  1.0,   },  //    x + y + z + w <= 20 .. b1
    { 10.0 , -2.0 , -1.0 , 1.0 ,  11.0,   },  //  -2x - y + z + w <= 10 .. b2
    { 10.0 ,  4.0 ,  1.0 , 0.0 , -1.0,   },  //    4x + y     - w <= 10 .. b3
  }
};


/*
lpp problem  = { 4, 4, {              // problem size 4*4
    {  0.0 ,  5.0 , -2.0 , 3.0 ,   },  // Max: z = 6*x + 4*y + 3*z,
    {  2.0 , -2.0 , -2.0 , 1.0 ,   },  //    -4*x + -5*y + -3z <= -40 .. b1
    {  3.0 ,  3.0 , -4.0 , 0.0 ,   },  //  -2x - 2y -6 z  <= -50 .. b2
    {  5.0 ,  0.0 ,  1.0 , 3.0 ,   },  //    -3x -4y -2z <= -60 .. b3
  }
};
*/

int main()
{
    print_problem_step(&problem,"Initial");
    solve_by_simplex(&problem);
    return 0;
}
