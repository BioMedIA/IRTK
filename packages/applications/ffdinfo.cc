/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkTransformation.h>

void usage()
{
    cerr << "Usage: ffdinfo [dofin]\n" << endl;
    cerr << "<-output file>       ffd value output file, prefix for each level, " << endl;
    cerr << "for example info (output will be info0.txt...infon.txt)" << endl;
    cerr << "<-sparsityoutput file>       sparsity value output file " << endl;
    exit(1);
}

char *resultout_name = NULL;
char *sparsityresultout_name = NULL;

int main(int argc, char **argv)
{
    int i, j, k, n, active, passive, ok;
    double dx, dy, dz, xaxis[3], yaxis[3], zaxis[3],l1_parameter,l1_deformation;

    // Check command line
    if (argc < 2) {
        usage();
    }

    l1_deformation = 0;
    l1_parameter = 0;

    // Read transformation
    irtkTransformation *transform = irtkTransformation::New(argv[1]);
    argv++;
    argc--;

    // Parse remaining parameters
    while (argc > 1){
        ok = false;
        if ((ok == false) && (strcmp(argv[1], "-output") == 0)){
            argc--;
            argv++;
            resultout_name = argv[1];
            argc--;
            argv++;
            ok = true;
        }
        if ((ok == false) && (strcmp(argv[1], "-sparsityoutput") == 0)){
            argc--;
            argv++;
            sparsityresultout_name = argv[1];
            argc--;
            argv++;
            ok = true;
        }
        if (ok == false){
            cerr << "Can not parse argument " << argv[1] << endl;
            usage();
        }
    } 

    irtkMultiLevelFreeFormTransformation *mffd = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(transform);
    cout << "Done" << endl;

    // Print transformation
    cout << "Transformation is " << mffd->NameOfClass() << endl;

    for (n = 0; n < mffd->NumberOfLevels(); n++) {

        // Extract current transformation level
        irtkFreeFormTransformation3D *ffd = dynamic_cast<irtkFreeFormTransformation3D *>(mffd->GetLocalTransformation(n));

        if (ffd == NULL) {
            cerr << "Free-form transformation is not 3D" << endl;
            exit(1);
        }

        ffd->GetOrientation(xaxis, yaxis, zaxis);
        ffd->GetSpacing(dx, dy, dz);

        // Print information about local transformation
        cout << "Local transformation no. " << n+1 << ":" << endl;
        cout << "Local transformation is a " << ffd->NameOfClass() << endl;
        cout << "Control points: \t" << ffd->GetX() << " x " << ffd->GetY()
            << " x " << ffd->GetZ() << endl;
        cout << "Orientation: \t\t" << xaxis[0] << " " << xaxis[1] << " "
            << xaxis[2] << endl;
        cout << "\t\t\t"  <<yaxis[0] << " " << yaxis[1] << " " << yaxis[2] << endl;
        cout << "\t\t\t" <<zaxis[0] << " " << zaxis[1] << " " << zaxis[2] << endl;
        cout << "Spacing: \t\t" << dx << " " << dy << " " << dx << endl;

        active  = 0;
        passive = 0;
        for (i = 0; i < ffd->GetX(); i++) {
            for (j = 0; j < ffd->GetY(); j++) {
                for (k = 0; k < ffd->GetZ(); k++) {
                    _Status sx, sy, sz;
                    ffd->GetStatus(i, j, k, sx, sy, sz);
                    if (sx == _Active) {
                        active++;
                    } else {
                        passive++;
                    }
                    if (sy == _Active) {
                        active++;
                    } else {
                        passive++;
                    }
                    if (sz == _Active) {
                        active++;
                    } else {
                        passive++;
                    }

                    dx = i;
                    dy = j;
                    dz = k;
                    if(sparsityresultout_name){
                        if(n == mffd->NumberOfLevels() - 1){
                            ffd->LatticeToWorld(dx,dy,dz);
                            mffd->LocalDisplacement(dx,dy,dz);
                            l1_deformation += fabs(dx)+fabs(dy)+fabs(dz);
                        }

                    }

                }
            }
        }
        cout << "Active control points:  " << active << " ("
            << active*100.0/(active+passive) << "%)" << endl;
        cout << "Passive control points: " << passive << " ("
            << passive*100.0/(active+passive) << "%)" << endl;

        if(resultout_name){
            char buffer[255];
            sprintf(buffer, "%s%d.txt",resultout_name,n);
            cerr << "Writing Results: " << buffer << endl;
            ofstream fout(buffer,ios::app);
            for(i = 0; i < ffd->NumberOfDOFs(); i++){
                fout << ffd->Get(i) << " ";
            }
            fout << endl;
            fout.close();
        }

        if(sparsityresultout_name){
            for(i = 0; i < ffd->NumberOfDOFs(); i++){
                l1_parameter +=  fabs(ffd->Get(i));
            }
        }
    }

    if(sparsityresultout_name){
        char buffer[255];
        sprintf(buffer, "%s",sparsityresultout_name);
        cerr << "Writing Sparsity Results: " << buffer << endl;
        ofstream fout(buffer,ios::app);  
        fout << l1_parameter << " " << l1_deformation << endl;
        fout.close();
    }
}
