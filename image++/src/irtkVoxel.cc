/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkVoxel.h>

std::string DataTypeName(int dataType)
{
	std::string typeName;

	switch (dataType){
	case IRTK_VOXEL_CHAR:
		typeName = "char";
		break;
	case IRTK_VOXEL_UNSIGNED_CHAR:
		typeName = "unsigned char";
		break;
	case IRTK_VOXEL_SHORT:
		typeName = "short";
		break;
	case IRTK_VOXEL_UNSIGNED_SHORT:
		typeName = "unsigned short";
		break;
	case IRTK_VOXEL_INT:
		typeName = "int";
		break;
	case IRTK_VOXEL_UNSIGNED_INT:
		typeName = "unsigned int";
		break;
	case IRTK_VOXEL_FLOAT:
		typeName = "float";
		break;
	case IRTK_VOXEL_DOUBLE:
		typeName = "double";
		break;
	case IRTK_VOXEL_RGB:
		typeName = "RGB";
		break;
	default:
		typeName = "unknown";
	}

	return typeName;
}
