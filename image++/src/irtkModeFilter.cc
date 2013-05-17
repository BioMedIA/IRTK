
#include <irtkImage.h>
#include <nr.h>

#include <irtkModeFilter.h>

#include <map>

#ifdef WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif



template <class VoxelType> irtkModeFilter<VoxelType>::irtkModeFilter()
{
	// Default connectivity.
	this->_Connectivity = CONNECTIVITY_26;
}

template <class VoxelType> irtkModeFilter<VoxelType>::~irtkModeFilter(void)
{
}

template <class VoxelType> bool irtkModeFilter<VoxelType>::RequiresBuffering(void)
{
  return true;
}

template <class VoxelType> const char *irtkModeFilter<VoxelType>::NameOfClass()
{
  return "irtkModeFilter";
}

template <class VoxelType> void irtkModeFilter<VoxelType>::Initialize()
{
  // Do the initial set up
  this->irtkImageToImage<VoxelType>::Initialize();

  this->_offsets.Initialize(this->_input, this->_Connectivity);
}



template <class VoxelType> void irtkModeFilter<VoxelType>::Run()
{
  int i, x, y, z, t, maskSize;
  VoxelType value;
  VoxelType *ptr2current, *ptr2offset;

  map<short, short> labelCount;
  map<short, short>::iterator iter;
  int ties, maxCount;
  short mode = 0;
  double inputMin, inputMax;

  int *tiedLabels;

  int randChoice;
  long ran2Seed;
  long ran2initialSeed;

  // Do the initial set up
  this->Initialize();


  // Get ready for random stuff.
  time_t tv;
  tv = time(NULL);

  ran2Seed = tv;
  ran2initialSeed = -1 * ran2Seed;
  (void) ran2(&ran2initialSeed);


  this->_input->GetMinMaxAsDouble(&inputMin, &inputMax);
  tiedLabels = new int[(int) inputMax];

  maskSize = this->_offsets.GetSize();

  for (t = 0; t < this->_input->GetT(); t++) {

    for (z = 0; z < this->_input->GetZ(); z++) {
      for (y = 0; y < this->_input->GetY(); y++) {
        this->_output->Put(0, y, z, t, this->_input->Get(0, y, z, t));
      }
    }

    for (z = 0; z < this->_input->GetZ(); z++) {
      for (x = 0; x < this->_input->GetX(); x++) {
        this->_output->Put(x, 0, z, t, this->_input->Get(x, 0, z, t));
      }
    }

    for (y = 0; y < this->_input->GetY(); y++) {
      for (x = 0; x < this->_input->GetX(); x++) {
        this->_output->Put(x, y, 0, t, this->_input->Get(x, y, 0, t));
      }
    }

    for (z = 1; z < this->_input->GetZ() - 1; z++) {
      for (y = 1; y < this->_input->GetY() - 1; y++) {
        for (x = 1; x < this->_input->GetX() - 1; x++) {

          labelCount.clear();

          ptr2current = this->_input->GetPointerToVoxels(x, y, z, t);

          // Collect labels from neighbourhood.
          for (i = 0; i < maskSize; ++i) {
            ptr2offset = ptr2current + this->_offsets(i);
            labelCount[*ptr2offset]++;
          }

          // Seek modal label (but there may be ties for the mode)
          maxCount = 0;
          for (iter = labelCount.begin(); iter != labelCount.end(); ++iter){
            if (iter->second > maxCount){
              maxCount = iter->second;
              mode     = iter->first;
            }
          }

          ties = 0;
          for (iter = labelCount.begin(); iter != labelCount.end(); ++iter){
            if (iter->second == maxCount){
              tiedLabels[ties] = iter->first;
              ++ties;
            }
          }

          if (ties > 1){
            randChoice = (int) floor( ran2(&ran2Seed) * ties );
            value = tiedLabels[ randChoice ];
          } else {
            value = mode;
          }

          this->_output->Put(x, y, z, t, value);
        }
      }
    }
  }

  delete [] tiedLabels;

  // Do the final cleaning up
  this->Finalize();
}

template class irtkModeFilter<irtkBytePixel>;
template class irtkModeFilter<irtkGreyPixel>;
