/* 
* MCCCS (Monte Carlo for Complex Chemical Systems) Towhee version*
* is a Gibbs ensemble Monte Carlo simulation code which is fully *
* functional for linear and branched molecules                   *
* with bond stretching, bond bending, dihedral potentials,       *
* 1-4 interactions, and intramolecular cross terms               *
* Copyright (C) 2004 Marcus G. Martin                            *
* see the file license.gpl for the full license information      *
*                                                                *
* This program is free software; you can redistribute it and/or  *
* modify it under the terms of the GNU General Public License    *
* as published by the Free Software Foundation; either version 2 *
* of the License, or (at your option) any later version.         *
*                                                                *
* This program is distributed in the hope that it will be useful,*
* but WITHOUT ANY WARRANTY; without even the implied warranty of *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
* GNU General Public License for more details.                   *
*                                                                *
* You should have received a copy of the GNU General Public      *
* License along with this program; if not, write to the Free     *
* Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,*
* MA  02111-1307, USA.                                           */

/*
  car2towhee takes a car format file and converts it into information
  needed to run Towhee

  Originally written 09-14-2004 by M.G. Martin
  Last modified 09-14-2007 by M.G. Martin
*/

/* 
   generic c includes 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* function declarations */
void clearline (FILE *);

int main()
{
  char filename[50];
  char atype[6],junk[6],prev[6]="go",current[6]="go",end[6]="end";
  int atomcount;
  float xxx,yyy,zzz,qqq;
  /* fp: input file descriptor */
  FILE *fp=NULL,*ofp=NULL;

  /* get the input file name */
  printf("Please enter a car filename to convert.\n");

  scanf("%s",filename);
  printf("Input filename is %s\n",filename);
  fp=fopen(filename,"r");
  printf("Please enter a car filename to convert.\n");

  /* set the output file name */
  strcpy(filename,"towhee_coords");
  ofp=fopen(filename,"w");

  /* skip the header information */
  clearline(fp);
  clearline(fp);
  clearline(fp);
  clearline(fp);
  /* get the box dimensions */
  fscanf(fp,"%s",junk);
  fscanf(fp,"%f",&xxx);
  fscanf(fp,"%f",&yyy);
  fscanf(fp,"%f",&zzz);
  clearline(fp);
  printf("Box dimensions %lf %lf %lf \n",xxx,yyy,zzz);

  /* initialize atom count */
  atomcount = 0;
  /* loop through the file */
  fscanf(fp,"%s",current);
  while ( strcmp(prev,end) != 0 || strcmp(current,end) != 0 ) {
    if ( strcmp(current,end) == 0 ) {
      /* just clear to the next line */
      clearline(fp);
    }
    else {
      /* read the coordinates and the atom name */
      fscanf(fp,"%f",&xxx);
      fscanf(fp,"%f",&yyy);
      fscanf(fp,"%f",&zzz);
      fscanf(fp,"%s",junk);
      fscanf(fp,"%s",junk);
      fscanf(fp,"%s",atype);
      fscanf(fp,"%s",junk);
      fscanf(fp,"%f",&qqq);
      clearline(fp);
      /* increment atomcount */
      atomcount++;
      /* write the coordinates and the atomname */
      fprintf(ofp,"%f %f %f %s %f \n",xxx,yyy,zzz,atype,qqq);
    }
    /* set prev to current for next loop */
    strcpy(prev,current);
    /* read in the next current */
    fscanf(fp,"%s",current);
  }
  /* close the files */
  fclose(fp);
  fclose(ofp);
  printf("towhee_coords successfully created \n");
  printf("Total number of atoms %d\n",atomcount);
}

/* 
   clearline is a simple utility to read to a newline or end of file and then return
   originally written 06-20-2003 by M.G. Martin
   last modified 06-20-2003 by M.G. Martin
*/
void clearline(FILE * fp)
{   
  int c;
  while ( ((c=getc(fp)) != EOF) && (c !='\n' ) &&  (c !='\r') );
}
