/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkImageHistograme_1D.cc 2 2008-12-23 12:40:14Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2008-12-23 12:40:14 +0000 (‰∫? 23 ÂçÅ‰∫åÊú?2008) $
  Version   : $Revision: 2 $
  Changes   : $Author: dr $

=========================================================================*/

#include <irtkImage.h>

#include <irtkHistogram.h>

template <class VoxelType> void irtkImageHistogram_1D<VoxelType>::Evaluate(irtkGenericImage<VoxelType> *image, double padding)
{  
	double value,min,max;
	int i,j,k,l;
	image->GetMinMaxAsDouble(&min,&max);
	this->PutMin(round(min));
	this->PutMax(round(max));
	this->PutNumberOfBins(round(max)-round(min)+1);
	for (l = 0; l < image->GetT(); l++){
		for (k = 0; k < image->GetZ(); k++){
			for (j = 0; j < image->GetY(); j++){
				for (i = 0; i < image->GetX(); i++){
					value = image->GetAsDouble(i, j, k, l);
					if(value > padding)
						this->AddSample(value);
				}
			}
		}
	}
}

template <class VoxelType> void irtkImageHistogram_1D<VoxelType>::BackProject(irtkGenericImage<VoxelType> *image)
{  
	VoxelType value;
	int i,j,k,l;
	for (l = 0; l < image->GetT(); l++){
		for (k = 0; k < image->GetZ(); k++){
			for (j = 0; j < image->GetY(); j++){
				for (i = 0; i < image->GetX(); i++){
					value = image->GetAsDouble(i, j, k, l);
					value = this->ValToBin(value);
					value = this->_bins[round(value)];
					image->PutAsDouble(i,j,k,l,value);
				}
			}
		}
	}
}

template <class VoxelType> void irtkImageHistogram_1D<VoxelType>::Equalize(VoxelType min,VoxelType max)
{
	int i;
	double count = 0;
	for(i=0;i<this->_nbins;i++){
		count += this->BinToPDF(i);
		this->_bins[i] = count*(max - min) + min;
	}
}

template class irtkImageHistogram_1D<unsigned char>;
template class irtkImageHistogram_1D<short>;
template class irtkImageHistogram_1D<unsigned short>;
template class irtkImageHistogram_1D<float>;
template class irtkImageHistogram_1D<double>;
