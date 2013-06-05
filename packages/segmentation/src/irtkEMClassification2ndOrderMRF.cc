/*
 * irtkEMClassification2ndOrderMRF.cc
 *
 *  Created on: 3 Jun 2013
 *      Author: cl6311
 */

/*
 * 	This class implements the Paper:
 * 	Neonatal brain segmentation using second order neighboorhood information, Ledig C. et al
 * 	PaPI Workshop 2012 in conjunction with MICCAI 2012, Nice *
 *
 */


#include <irtkImage.h>
#include <vector>
#include <utility>

#include <set>

#include <irtkEMClassification2ndOrderMRF.h>
#include <irtkGaussianBlurring.h>


// Default constructor
irtkEMClassification2ndOrderMRF::irtkEMClassification2ndOrderMRF() : irtkEMClassification()
{
  cerr<<"irtkEMClassification2ndOrderMRF() ";
}

irtkEMClassification2ndOrderMRF::irtkEMClassification2ndOrderMRF(int noTissues, irtkRealImage **atlas, irtkRealImage *background) : irtkEMClassification(noTissues, atlas, background)
{
	_isLogTransformed = false;
	_has_background = true;
	_has_MRF_2nd_order = false;
	SetMRFWeights(0.15, 0.75);
	SetRelaxationFactor(0.5);
	SetRelaxationFactor_2nd(0.5);
}

irtkEMClassification2ndOrderMRF::irtkEMClassification2ndOrderMRF(int noTissues, irtkRealImage **atlas ) : irtkEMClassification(noTissues, atlas )
{
	_isLogTransformed = false;
	_has_background = false;
	_has_MRF_2nd_order = false;
	SetMRFWeights(0.15, 0.75);
	SetRelaxationFactor(0.5);
	SetRelaxationFactor_2nd(0.5);
}

void irtkEMClassification2ndOrderMRF::SetMRFWeights( double adjacent, double distant )
{
	_mrf_weight_adjacent = adjacent;
	_mrf_weight_distant = distant;
}

void irtkEMClassification2ndOrderMRF::SetRelaxationFactor( double relax_factor )
{
	_relax_factor = relax_factor;
}

void irtkEMClassification2ndOrderMRF::SetRelaxationFactor_2nd( double relax_factor )
{
	_relax_factor_2nd = relax_factor;
}

void irtkEMClassification2ndOrderMRF::SetMRF_2nd_order(vector< vector< pair< pair<int,int>, double > > > &mrf_2nd)
{
	_has_MRF_2nd_order = true;
	_connectivity_2nd_order = mrf_2nd;
}

void irtkEMClassification2ndOrderMRF::SetLogTransformed(bool state)
{
	_isLogTransformed = state;
}


void irtkEMClassification2ndOrderMRF::SetInput(const irtkRealImage &image, const irtkMatrix &connectivity)
{
  _uncorrected = image;
  if( connectivity.Cols() != _number_of_tissues || connectivity.Rows() != connectivity.Cols() )
  {
	  cerr << "Warning: Connectivity matrix has wrong size! expected: " << _number_of_tissues << "x"<< _number_of_tissues << endl;
	  cerr << "is:" << connectivity.Rows() << "x" << connectivity.Rows() << endl;
  }
  _connectivity = connectivity;

  irtkEMClassification::SetInput(_uncorrected);


  // Initialize MRF field
  _MRF_weights = _uncorrected;
  irtkRealPixel *ptr = _MRF_weights.GetPointerToVoxels();
  irtkRealPixel *pm = _mask.GetPointerToVoxels();

  for( int i = 0; i < _uncorrected.GetNumberOfVoxels(); ++i )
  {
	  if( *pm == 1 )
	  {
		  *ptr = 1.0;
	  }
	  pm++;
	  ptr++;
  }
}


void irtkEMClassification2ndOrderMRF::BStep()
{
  // Create bias correction filter
  _biascorrection.SetInput(&_uncorrected, &_estimate);
  _biascorrection.SetWeights(&_weights);

  _biascorrection.SetOutput(_biasfield);
  _biascorrection.SetPadding((short int) _padding);
  _biascorrection.SetMask(&_mask);
  _biascorrection.Run();

  // Generate bias corrected image for next iteration
  _input = _uncorrected;
  _biascorrection.Apply(_input);
}

void irtkEMClassification2ndOrderMRF::removePVclasses(double threshold )
{
	irtkRealPixel *pm = _mask.GetPointerToVoxels();
	irtkRealPixel *ptr = _input.GetPointerToVoxels();
	_output.First();

	double* tempResult = new double[_number_of_tissues];

	for( int i = 0; i < _input.GetNumberOfVoxels(); ++i )
	{
		if( *pm == 1 )
		{
			for( int k = 0; k < _number_of_tissues; ++k )
			{
				tempResult[k] = (_output.GetValue(k) > threshold ) ? _output.GetValue(k) : 0;
			}
			map<int,int>::iterator iter = pv_classes.begin();
			while( iter != pv_classes.end() )
			{
				int position = iter->second;
				int classA = pv_connections[position].first;
				int classB = pv_connections[position].second;
				int classPV = iter->first;
				double miA = _mi[classA];
				double miB = _mi[classB];

				if( _output.GetValue(classPV) > threshold )
				{
					if( _output.GetValue(classA) < _output.GetValue(classPV) && _output.GetValue(classB) < _output.GetValue(classPV) )
					{
						double expec_sum = _output.GetValue(classA) + _output.GetValue(classB) + _output.GetValue(classPV);
						double frac_content=(miA-(*ptr))/(miA-miB);
						frac_content = (frac_content < 0) ? 0 : frac_content;
						frac_content = (frac_content > 1) ? 1 : frac_content;
						tempResult[classA] = (1.0-frac_content) / expec_sum;
						tempResult[classB] = frac_content / expec_sum;
					}
					else if( _output.GetValue(classA) > _output.GetValue(classB) )
					{
						tempResult[classA] = 1;
						tempResult[classB] = 0;
					}
					else
					{
						tempResult[classA] = 0;
						tempResult[classB] = 1;
					}
					_output.SetValue(classPV,0);
				}
				tempResult[classPV] = 0;
				iter++;
			}
			double denom = 0;
			for( int k = 0; k < _number_of_tissues; ++k )
			{
				denom += tempResult[k];
			}
			for( int k = 0; k < _number_of_tissues; ++k )
			{
				_output.SetValue(k, tempResult[k] / denom );
			}
		}
		_output.Next();
		ptr++;
		pm++;
	}
	delete[] tempResult;
}

void irtkEMClassification2ndOrderMRF::RStep_2nd_order()
{
	cout << "Prior relaxation based on 2nd order connectivity" << endl;
	double relaxFactor = _relax_factor_2nd;

	cout << endl << relaxFactor << endl;
	double dx,dy,dz;
	_input.GetPixelSize(&dx,&dy,&dz);

	double sx, sy,sz, st;
	sx = 1.0/dx;
	sy = 1.0/dy;
	sz = 1.0/dz;
	// do not do spatial smoothing !!!
	st = 0.0;

	double* numerator = new double[_number_of_tissues];
	double* cumulated_probs = new double[_number_of_tissues];

	irtkRealPixel *pm = _mask.GetPointerToVoxels();
	int per = 0;

	for (int i=0; i< _input.GetNumberOfVoxels(); i++)
	{
	    if (i*10.0/_input.GetNumberOfVoxels() > per) {
	      per++;
	      cerr<<per<<"0%...";
	    }
	    int index = i;
		double denominator = 0.0;
		if (*pm == 1)
		{

			int x,y,z,t;

			t = index / (_input.GetX() * _input.GetY() * _input.GetZ());
			index -= t * (_input.GetX() * _input.GetY() * _input.GetZ());

		    z = index / (_input.GetX() * _input.GetY());
		    index -= z * (_input.GetX() * _input.GetY());

		    x = index % _input.GetX();
		    y = index / _input.GetX();

			int lt = max(t-1,0);
			int rt = min(t+1, _input.GetT()-1 );
			int lx = max(x-1,0);
			int rx = min(x+1, _input.GetX()-1 );
			int ly = max(y-1,0);
			int ry = min(y+1, _input.GetY()-1 );
			int lz = max(z-1,0);
			int rz = min(z+1, _input.GetZ()-1 );

			for( int tissue = 0; tissue < _number_of_tissues; ++tissue )
			{
				vector< pair< pair<int,int>, double > > vec = _connectivity_2nd_order[tissue];
				vector< pair< pair<int,int>, double > >::iterator iter = vec.begin();
				if( iter == vec.end() )
				{
					continue;
				}

				for( int k = 0; k < _number_of_tissues; k++)
				{
					cumulated_probs[k] = 0;
					cumulated_probs[k] += sx*(_output.GetValue( rx ,y,z,t,k) + _output.GetValue( lx,y,z,t,k));
					cumulated_probs[k] += sy*(_output.GetValue(x,ry,z,t,k) + _output.GetValue(x,ly,z,t,k));
					cumulated_probs[k] += sz*(_output.GetValue(x,y,rz,t,k) + _output.GetValue(x,y,lz,t,k));

					if( rt != t )
						cumulated_probs[k] += st*_output.GetValue(x,y,z,rt,k);
					if( lt != t )
						cumulated_probs[k] += st*_output.GetValue(x,y,z,lt,k);
				}

				for( ; iter != vec.end(); ++iter )
				{
					pair<int,int> neighbors = (*iter).first;
					int neighborA = neighbors.first;
					int neighborB = neighbors.second;

					double intersect = min(cumulated_probs[neighborA],cumulated_probs[neighborB]);

					if( intersect > 0.5 )
					{
						double tmp1 = _atlas.GetValue(x,y,z,t,tissue);
						double tmp2 = _atlas.GetValue(x,y,z,t,neighborA);
						double tmp3 = _atlas.GetValue(x,y,z,t,neighborB);
						_atlas.SetValue(x,y,z,t,tissue,tmp1 * (1.0 - relaxFactor));
						_atlas.SetValue(x,y,z,t,neighborA,tmp2+tmp1*relaxFactor/2.0);
						_atlas.SetValue(x,y,z,t,neighborB,tmp3+tmp1*relaxFactor/2.0);
					}
				}
			}
			for( int k = 0; k < _number_of_tissues; ++k )
			{
				numerator[k] = _atlas.GetValue(x,y,z,t,k);
				denominator += numerator[k];
			}

			for( int k = 0; k < _number_of_tissues; ++k )
			{
				if( denominator != 0 )
				{
					_atlas.SetValue(x,y,z,t,k, numerator[k] / denominator);
				}
				else
				{
					cerr << "Division by 0 while computing relaxed prior probabilities" << endl;
					if (k==0) _atlas.SetValue(x,y,z,t,k, 1);
					else     _atlas.SetValue(x,y,z,t,k, 0);
				}
			}
		}
		pm++;
	}
	delete[] numerator;
	delete[] cumulated_probs;
}

// Relaxation according to Cardoso in MICCAI 2011
void irtkEMClassification2ndOrderMRF::RStep()
{
	double relaxFactor = _relax_factor;

	_atlas.First();
	_output.First();
	irtkProbabilisticAtlas filteredAtlas;

	irtkRealImage filteredImage;
	irtkRealImage filterInput;
	cerr << "Filtering: ";
	for( int k = 0; k < _number_of_tissues; ++k )
	{
		cerr << k << " ";
		filteredImage = _output.GetImage(k);
		filterInput = _output.GetImage(k);
		irtkGaussianBlurring<irtkRealPixel> filter(2.0);
		filter.SetInput(&filterInput);
		filter.SetOutput(&filteredImage);
		filter.Run();
		filteredAtlas.AddImage(filteredImage);
	}
	cerr << endl;

	filteredAtlas.First();
	irtkRealPixel *ptr = _input.GetPointerToVoxels();
	irtkRealPixel *pm = _mask.GetPointerToVoxels();

	int per = 0;
	double* numerator = new double[_number_of_tissues];
	for (int i=0; i< _input.GetNumberOfVoxels(); i++)
	{
	    if (i*10.0/_input.GetNumberOfVoxels() > per) {
	      per++;
	      cerr<<per<<"0%...";
	    }
		double denominator = 0.0;
		if (*pm == 1)
		{
			for( int k = 0; k < _number_of_tissues; ++k )
			{
				filteredAtlas.SetValue(k,  relaxFactor * filteredAtlas.GetValue(k) + (1.0-relaxFactor) * _atlas.GetValue(k) );
			}

			for( int k = 0; k < _number_of_tissues; ++k )
			{
				numerator[k] = .0;
				double temp = .0;
				for( int j = 0; j < _number_of_tissues; ++j )
				{
					if ( _connectivity(k,j) <= 1 )
					{
						temp += filteredAtlas.GetValue(k) * filteredAtlas.GetValue(j);
					}
				}
				numerator[k] += temp;
				denominator += numerator[k];
			}

			for( int k = 0; k < _number_of_tissues; ++k )
			{
				if( denominator != 0 )
				{
					_atlas.SetValue(k, numerator[k] / denominator);
				}
				else
				{
					cerr << "Division by 0 while computing relaxed prior probabilities" << endl;
					if (k==0) _atlas.SetValue(k, 1);
					else     _atlas.SetValue(k, 0);
				}
			}
		}
		pm++;
		ptr++;
		_output.Next();
		_atlas.Next();
		filteredAtlas.Next();
	}
	delete[] numerator;
}


int irtkEMClassification2ndOrderMRF::AddPartialVolumeClass(int classA, int classB)
{
	irtkRealImage pvclass = _input;
	irtkRealPixel *ptr_pvclass = pvclass.GetPointerToVoxels();

	// mixing coefficient
	double gamma = 0.0;
	int N = 0;
	irtkRealPixel *ptr = _input.GetPointerToVoxels();
	irtkRealPixel *pm = _mask.GetPointerToVoxels();

	for( int i = 0; i < pvclass.GetNumberOfVoxels(); ++i )
	{
		if( *pm == 1 )
		{
			double fc = (_mi[classA] - *ptr ) / ( _mi[classA]-_mi[classB]);
			if( fc >= 0 && fc <= 1.0 )
			{
				gamma += fc;
				N++;
			}
		}
		*ptr_pvclass = .0;
		ptr_pvclass++;
		pm++;
		ptr++;
	}
	if( N > 0 )
	{
		gamma /= N;
	}
	else
	{
		cerr << "No mixel voxels found, not adding partial volume class!" << endl;
		return -1;
	}

	double* mi = new double[_number_of_tissues+1];
	double* sigma = new double[_number_of_tissues+1];

	mi[_number_of_tissues] = 0;
	sigma[_number_of_tissues] = 0;
	for( int i = 0; i < _number_of_tissues; ++i )
	{
		mi[i] = _mi[i];
		sigma[i] = _sigma[i];
	}
	mi[_number_of_tissues] = (1.0 - gamma) * mi[classA] + gamma * mi[classB];
	sigma[_number_of_tissues] = (1.0 - gamma) * (1.0 - gamma) * sigma[classA] + gamma * gamma * sigma[classB];

	_atlas.First();
	_output.First();
	ptr = pvclass.GetPointerToVoxels();
	pm = _mask.GetPointerToVoxels();

	for( int i = 0; i < _input.GetNumberOfVoxels(); ++i )
	{
		if( *pm == 1 )
		{
			double tmp = _output.GetValue(classA) * _output.GetValue(classB);
			if( tmp > 0.0 )
			{
				tmp = sqrt(tmp) / 0.5;
			}
			else
			{
				tmp = 0.0;
			}
			*ptr = tmp;
			for( int k = 0; k < _number_of_tissues; ++k )
			{
				tmp += _output.GetValue(k);
			}
			if( tmp > 0.0 )
			{
				*ptr /= tmp;
				for( int k = 0; k < _number_of_tissues; ++k )
				{
					_atlas.SetValue(k,_output.GetValue(k) / tmp);
				}
			}
			else
			{
				cerr << "error probability = 0" << endl;
				return -1;
			}

		}
		ptr++;
		pm++;
		_atlas.Next();
		_output.Next();
	}

	delete[] _mi;
	delete[] _sigma;

	_number_of_tissues++;
	_mi = new double[_number_of_tissues];
	_sigma = new double[_number_of_tissues];
	for( int k = 0; k < _number_of_tissues; ++k)
	{
		_mi[k] = mi[k];
		_sigma[k] = sigma[k];
	}
	delete[] mi;
	delete[] sigma;

	_atlas.AddImage(pvclass);

	irtkRealImage newimage = pvclass;
	_output.AddImage(newimage);

	_atlas.First();
	_output.First();
	for( int i = 0; i < _input.GetNumberOfVoxels(); ++i )
	{
		for( int k = 0; k < _number_of_tissues; ++k)
		{
			_output.SetValue(k, _atlas.GetValue(k));
		}
		_atlas.Next();
		_output.Next();
	}
	cout << "Connectivity before update" << endl;
	_connectivity.Print();
	// PV classes get same connectivity as parent class!
	irtkMatrix newconnectivity(_connectivity.Rows()+1, _connectivity.Cols()+1);
	for( int i = 0; i < newconnectivity.Rows(); ++i )
	{
		for( int j = 0; j < newconnectivity.Cols(); ++j )
		{
			if( i < _connectivity.Rows() && j < _connectivity.Cols() )
			{
				// connectivity basically stays at it was!
				newconnectivity.Put(i,j,_connectivity.Get(i,j));
			}
			// same class
			else if( i == j )
			{
				newconnectivity.Put(i,j,0);
			}
			// this is for the PV class: close to pv contributing classes, distant to all others
			else
			{
				if( i == classA || i == classB || j == classA || j == classB )
				{
					newconnectivity.Put(i,j, 1);
				}
				else
				{
					newconnectivity.Put(i,j, 2);
				}
			}
		}
	}

	// if classification is done with background ensure that background class is the last tissue class
	// -> swap prob maps and adapt connectivity matrix
	if( _has_background )
	{
		double tmp =_mi[_number_of_tissues-1];
		_mi[_number_of_tissues-1] = _mi[_number_of_tissues-2];
		_mi[_number_of_tissues-2] = tmp;
		_sigma[_number_of_tissues-1] = _sigma[_number_of_tissues-2];
		_sigma[_number_of_tissues-2] = tmp;

		_output.SwapImages(_number_of_tissues-1, _number_of_tissues-2);
		_atlas.SwapImages(_number_of_tissues-1, _number_of_tissues-2);

		// now swap last two rows and columns of newconnectivity, since prob maps (for background) were swapped as well
		for( int i = 0; i < newconnectivity.Cols(); ++i )
		{
			double tmp = newconnectivity.Get(i,newconnectivity.Rows()-1);
			newconnectivity.Put(i,newconnectivity.Rows()-1, newconnectivity.Get(i,newconnectivity.Rows()-2));
			newconnectivity.Put(i,newconnectivity.Rows()-2, tmp);
		}
		for( int i = 0; i < newconnectivity.Rows(); ++i )
		{
			double tmp = newconnectivity.Get(newconnectivity.Cols()-1,i);
			newconnectivity.Put(newconnectivity.Cols()-1,i, newconnectivity.Get(newconnectivity.Cols()-2,i));
			newconnectivity.Put(newconnectivity.Cols()-2, i, tmp);
		}
	}

	_connectivity = newconnectivity;

	int pv_position = _number_of_tissues - 1;
	if( _has_background )
	{
		pv_position = _number_of_tissues - 2;
	}
	pv_classes.insert(make_pair(pv_position, pv_connections.size() ) );

	pv_connections.push_back( make_pair(classA, classB) );
	pv_fc.push_back(gamma);
	cout << "Fractional Content of classA=" << 1.0-gamma << endl;
	cout << "connectivity after update " << endl;
	_connectivity.Print();

	if( _has_MRF_2nd_order )
	{
		_connectivity_2nd_order.resize(_number_of_tissues);
		vector< pair< pair<int,int>, double> >::iterator iter;

		if( _has_background )
		{
			_connectivity_2nd_order[_number_of_tissues-1] = _connectivity_2nd_order[_number_of_tissues-2];
			for ( int k = 0; k < _number_of_tissues; ++k )
			{
				int conn_size = _connectivity_2nd_order[k].size();
				for( int l = 0; l < conn_size; ++l )
				{
					pair<int,int> neighbors = _connectivity_2nd_order[k][l].first;
					double w = _connectivity_2nd_order[k][l].second;
					int neighborA = neighbors.first;
					int neighborB = neighbors.second;
					_connectivity_2nd_order[k].erase(_connectivity_2nd_order[k].begin() + l);
					if( neighborA == _number_of_tissues-2 )
					{
						neighborA = _number_of_tissues-1;
					}
					if( neighborB == _number_of_tissues-2 )
					{
						neighborB = _number_of_tissues-1;
					}
					_connectivity_2nd_order[k].push_back( make_pair( make_pair(neighborA, neighborB), w ) );
				}
			}
		}
		// create connectivity for partial volume class as union of the participating classes
		_connectivity_2nd_order[pv_position].clear();
		for ( int k = 0; k < _number_of_tissues; ++k )
		{
			if( k == classA || k == classB )
			{
				int conn_size = _connectivity_2nd_order[k].size();
				for( int l = 0; l < conn_size; ++l )
				{
					pair<int,int> neighbors = _connectivity_2nd_order[k][l].first;
					double w = _connectivity_2nd_order[k][l].second;
					int neighborA = neighbors.first;
					int neighborB = neighbors.second;

					_connectivity_2nd_order[pv_position].push_back( make_pair( make_pair(neighborA, neighborB), w ) );
				}
			}
		}
		for ( int k = 0; k < _number_of_tissues; ++k )
		{
			int conn_size = _connectivity_2nd_order[k].size();
			for( int l = 0; l < conn_size; ++l )
			{
				pair<int,int> neighbors = _connectivity_2nd_order[k][l].first;
				double w = _connectivity_2nd_order[k][l].second;
				int neighborA = neighbors.first;
				int neighborB = neighbors.second;

				if( neighborA == classA ) neighborA = pv_position;
				if( neighborB == classB ) neighborB = pv_position;
				if( neighborA == pv_position || neighborB == pv_position )
				{
					_connectivity_2nd_order[k].push_back( make_pair( make_pair(neighborA, neighborB), w ) );
				}
			}
		}
	}
	return pv_position;
}

void irtkEMClassification2ndOrderMRF::RefineSegmentation(int bg = -1)
{
	  irtkRealImage segmentation;

	  int i, j;
	  double max_val;
	  irtkRealPixel *ptr;

	  cerr<<"Constructing segmentation"<<endl;

	  // Initialize pointers of probability maps
	  _output.First();

	  // Initialize segmentation to same size as input
	  segmentation = _input;
	  ptr = segmentation.GetPointerToVoxels();

	  irtkRealPixel *pm = _mask.GetPointerToVoxels();

	  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
	    if (*ptr != _padding && *pm == 1) {
	      max_val  = 0;
	      *ptr = 0;
	      for (j = 0; j < _number_of_tissues; j++) {
	        if (_output.GetValue(j) > max_val) {
	          max_val  = _output.GetValue(j);
	          *ptr = j;
	          //*ptr = j+1;
	          if ( (j+1) == _number_of_tissues && _has_background) *ptr=0;
	        }
	      }
	    }
	    else
	    {
	    	*ptr = 0;
	    }
	    pm++;
	    ptr++;
	    _output.Next();
	  }
	  _output.First();

	  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
			int x,y,z,t;
			int index = i;
			t = index / (_input.GetX() * _input.GetY() * _input.GetZ());
			index -= t * (_input.GetX() * _input.GetY() * _input.GetZ());

		    z = index / (_input.GetX() * _input.GetY());
		    index -= z * (_input.GetX() * _input.GetY());

		    x = index % _input.GetX();
		    y = index / _input.GetX();

			int lx = max(x-1,0);
			int rx = min(x+1, _input.GetX()-1 );
			int ly = max(y-1,0);
			int ry = min(y+1, _input.GetY()-1 );
			int lz = max(z-1,0);
			int rz = min(z+1, _input.GetZ()-1 );

			int label = segmentation.Get(x,y,z);
			int same_label = 0;

			int* cumulated_segs = new int[_number_of_tissues];
			memset(cumulated_segs, 0, sizeof(cumulated_segs));

			cumulated_segs[(int) segmentation.Get(lx,y,z)]++;
			cumulated_segs[(int) segmentation.Get(rx,y,z)]++;
			cumulated_segs[(int) segmentation.Get(x,ly,z)]++;
			cumulated_segs[(int) segmentation.Get(x,ry,z)]++;
			cumulated_segs[(int) segmentation.Get(x,y,lz)]++;
			cumulated_segs[(int) segmentation.Get(x,y,rz)]++;

			if( segmentation.Get(lx,y,z) == label ) same_label++;
			if( segmentation.Get(rx,y,z) == label ) same_label++;
			if( segmentation.Get(x,ly,z) == label ) same_label++;
			if( segmentation.Get(x,ry,z) == label ) same_label++;
			if( segmentation.Get(x,y,lz) == label ) same_label++;
			if( segmentation.Get(x,y,rz) == label ) same_label++;

			int max_label = 0;
			for ( int k = 0; k < _number_of_tissues; ++k )
			{
				if( cumulated_segs[max_label] < cumulated_segs[k] ) max_label = k;
			}

			bool refined = false;
			if( bg != -1 && label == bg && _output.GetValue(x,y,z,label) < 0.9 )
			{
				double bg_label = _output.GetValue(x,y,z,label);
				int nonzeros = 0;
				for( int k = 0; k < _number_of_tissues; ++k )
				{
					if( _output.GetValue(x,y,z,k) && k != bg ) nonzeros++;
				}
				for( int k = 0; k < _number_of_tissues; ++k )
				{
					double val_label = _output.GetValue(x,y,z,k);
					if( val_label )
					{
						_output.SetValue(x,y,z,k,val_label + bg_label / (nonzeros));
					}
				}
				_output.SetValue(x,y,z,label,0);
				refined = true;
			}

			if( same_label <= 0 && !isPVclass(label) && !refined && label != max_label )
			{
				double val_label = _output.GetValue(x,y,z,label);
				_output.SetValue(x,y,z,max_label, 0.5 * val_label + _output.GetValue(x,y,z,max_label));
				_output.SetValue(x,y,z,label, 0.5 * val_label );
			}
			delete[] cumulated_segs;
	  }
}

void irtkEMClassification2ndOrderMRF::ConstructSegmentation(irtkRealImage &segmentation)
{
  int i, j;
  double max;
  irtkRealPixel *ptr;

  cerr<<"Constructing segmentation"<<endl;

  // Initialize pointers of probability maps
  _output.First();

  // Initialize segmentation to same size as input
  segmentation = _input;
  ptr = segmentation.GetPointerToVoxels();

  irtkRealPixel *pm = _mask.GetPointerToVoxels();

  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
    if (*ptr != _padding && *pm == 1) {
      max  = 0;
      *ptr = 0;
      for (j = 0; j < _number_of_tissues; j++) {
        if (_output.GetValue(j) > max) {
          max  = _output.GetValue(j);
          *ptr = j;
          //*ptr = j+1;
          if ( (j+1) == _number_of_tissues && _has_background) *ptr=0;
        }
      }
    }
    else
    {
    	*ptr = 0;
    }
    pm++;
    ptr++;
    _output.Next();
  }
}

void irtkEMClassification2ndOrderMRF::ConstructSegmentationPV(irtkRealImage &segmentation)
{
  int i, j;
  double max;
  irtkRealPixel *pm = _mask.GetPointerToVoxels();
  irtkRealPixel *ptr;
  irtkRealPixel *ptrimage = _input.GetPointerToVoxels();

  cerr<<"Constructing segmentation"<<endl;

  // Initialize pointers of probability maps
  _output.First();

  // Initialize segmentation to same size as input
  segmentation = _input;
  ptr = segmentation.GetPointerToVoxels();

  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
	*ptr = 0;
	if (*ptr != _padding && *pm == 1) {
	  max  = 0;
	  for (j = 0; j < _number_of_tissues; j++) {
		double new_max = _output.GetValue(j);
		for( int k = j + 1; k < _number_of_tissues; ++k )
		{
			if( pv_classes.find(k) != pv_classes.end() )
			{
				map<int,int>::iterator iter =  pv_classes.find(k);
				int position = iter->second;
				double value = _output.GetValue(k);
				int classA = pv_connections[position].first;
				int classB = pv_connections[position].second;

				double avg_mi = (_mi[classA] + _mi[classB]) / 2.0;
				if( classA == j || classB == j )
				{
					if( (*ptrimage >= avg_mi && _mi[j] >= avg_mi) || (*ptrimage <= avg_mi && _mi[j] <= avg_mi) )
					{
						new_max += value;
					}
				}
			}
		}
		if ( new_max > max && pv_classes.find(j) == pv_classes.end()) {
		  max  = new_max;
		  *ptr = j+1;
		  if ( (j+1) == _number_of_tissues && _has_background) *ptr=0;
		}
	  }
	}
	pm++;
	ptrimage++;
	ptr++;
	_output.Next();
  }
}

double irtkEMClassification2ndOrderMRF::getMRFenergy_2nd_order(int index, int tissue)
{
	vector< pair< pair<int,int>, double > > vec = _connectivity_2nd_order[tissue];
	vector< pair< pair<int,int>, double > >::iterator iter = vec.begin();
	if( iter == vec.end() )
	{
		return 1.0;
	}

	double dx,dy,dz;
	_input.GetPixelSize(&dx,&dy,&dz);

	double sx, sy,sz, st;
	sx = 1.0/dx;
	sy = 1.0/dy;
	sz = 1.0/dz;
	// do not do spatial smoothing !!!
	st = 0.0;

	int x,y,z,t;

	t = index / (_input.GetX() * _input.GetY() * _input.GetZ());
	index -= t * (_input.GetX() * _input.GetY() * _input.GetZ());

    z = index / (_input.GetX() * _input.GetY());
    index -= z * (_input.GetX() * _input.GetY());

    x = index % _input.GetX();
    y = index / _input.GetX();

	int lt = max(t-1,0);
	int rt = min(t+1, _input.GetT()-1 );
	int lx = max(x-1,0);
	int rx = min(x+1, _input.GetX()-1 );
	int ly = max(y-1,0);
	int ry = min(y+1, _input.GetY()-1 );
	int lz = max(z-1,0);
	int rz = min(z+1, _input.GetZ()-1 );

	double weight_2nd = 0.0;

	double* cumulated_probs = new double[_number_of_tissues];

	for( int k = 0; k < _number_of_tissues; k++)
	{
		cumulated_probs[k] = 0;
		cumulated_probs[k] += sx*(_output.GetValue( rx ,y,z,t,k) + _output.GetValue( lx,y,z,t,k));
		cumulated_probs[k] += sy*(_output.GetValue(x,ry,z,t,k) + _output.GetValue(x,ly,z,t,k));
		cumulated_probs[k] += sz*(_output.GetValue(x,y,rz,t,k) + _output.GetValue(x,y,lz,t,k));

		if( rt != t )
			cumulated_probs[k] += st*_output.GetValue(x,y,z,rt,k);
		if( lt != t )
			cumulated_probs[k] += st*_output.GetValue(x,y,z,lt,k);
	}

	for( ; iter != vec.end(); ++iter )
	{
		pair<int,int> neighbors = (*iter).first;
		double w = (*iter).second;
		int neighborA = neighbors.first;
		int neighborB = neighbors.second;

		weight_2nd += w * cumulated_probs[neighborA] * cumulated_probs[neighborB];
	}

	delete[] cumulated_probs;

	double expo = -1.0 * weight_2nd;
	return exp(expo);
}

double irtkEMClassification2ndOrderMRF::getMRFenergy(int index, int tissue)
{
	if( _connectivity.Rows() == 1 )
	{
		return 1.0;
	}
	double dx,dy,dz;
	_input.GetPixelSize(&dx,&dy,&dz);

	double sx, sy,sz, st;
	sx = 1.0/dx;
	sy = 1.0/dy;
	sz = 1.0/dz;
	// do not do spatial smoothing !!!
	st = 0.0;

	double expo = -1.0;

	double energy = .0;
	int x,y,z,t;

	t = index / (_input.GetX() * _input.GetY() * _input.GetZ());
	index -= t * (_input.GetX() * _input.GetY() * _input.GetZ());

    z = index / (_input.GetX() * _input.GetY());
    index -= z * (_input.GetX() * _input.GetY());

    x = index % _input.GetX();
    y = index / _input.GetX();

	// if MRF field is spatially invariant keep this beta constant!
	// double beta = 1.0;
	double beta = _MRF_weights.Get(x,y,z);
	expo *= beta;

	// PAPER
//	double distant = 3.0;
//	double close = 0.5;
	// IMPLEMENTATION
	double distant = _mrf_weight_distant;
	double close = _mrf_weight_adjacent;

	double weight = 0;
	int lt = max(t-1,0);
	int rt = min(t+1, _input.GetT()-1 );
	int lx = max(x-1,0);
	int rx = min(x+1, _input.GetX()-1 );
	int ly = max(y-1,0);
	int ry = min(y+1, _input.GetY()-1 );
	int lz = max(z-1,0);
	int rz = min(z+1, _input.GetZ()-1 );

	for( int k = 0; k < _number_of_tissues; k++)
	{
		double temp = 0;
		weight = 0;
		double conn = _connectivity.Get(k, tissue);

		if( conn  == 1 )
			weight = close;
		else if( conn == 2 )
			weight = distant;

		temp += sx * (_output.GetValue( rx ,y,z,t,k) + _output.GetValue( lx,y,z,t,k) );
		temp += sy * ( _output.GetValue(x,ry,z,t,k) + _output.GetValue(x,ly,z,t,k) );
		temp += sz * ( _output.GetValue(x,y,rz,t,k) + _output.GetValue(x,y,lz,t,k) );
		if( lt != t )
			temp += st * _output.GetValue(x,y,z,rt,k);
		if( rt != t )
			temp += st * _output.GetValue(x,y,z,lt,k);
		energy += temp * weight;
	}

	expo *= energy;
	return exp(expo);
}

void irtkEMClassification2ndOrderMRF::EStepMRF_2nd_order()
{
	cout << "E-Step with 2nd order MRF" << endl;
	  int i, k;
	  double x;
	  irtkGaussian* G = new irtkGaussian[_number_of_tissues];

	  for (k = 0; k < _number_of_tissues; k++) {
	    G[k].Initialise( _mi[k], _sigma[k]);
	  }

	  _atlas.First();
	  _output.First();
	  irtkRealPixel *ptr = _input.GetPointerToVoxels();
	  irtkRealPixel *pm = _mask.GetPointerToVoxels();
	  int per = 0;
	  double* numerator = new double[_number_of_tissues];
	  double denominator=0, temp=0;
	  bool bMRF = _number_of_tissues == _connectivity.Rows();

	  if(!bMRF)
	  {
		  cout << "Warning: number of tissues does not match size of connectivity matrix!" << endl;
	  }
	  for (i=0; i< _input.GetNumberOfVoxels(); i++) {
		if (i*10.0/_input.GetNumberOfVoxels() > per) {
		  per++;
		  cerr<<per<<"0%...";
		}
		denominator = 0;
		temp = 0;

	    if (*pm == 1) {
	      x = *ptr;
	      double hp = 0.5;
		  double hm = 0.5;
		  if( _isLogTransformed )
		  {
				hp = log( exp(x)+0.5 ) - x;
				hm = x - log( (exp(x) - 0.5) > 0 ? exp(x) - 0.5 : exp(x) );
		  }
	      double* MRFenergies = new double[_number_of_tissues];

	      double denominatorMRF = .0;

	      for (k = 0; k < _number_of_tissues; k++) {
    		  MRFenergies[k] = _atlas.GetValue(k) * getMRFenergy(i,k) * getMRFenergy_2nd_order(i,k);
	    	  denominatorMRF += MRFenergies[k];
	      }

	      for (k = 0; k < _number_of_tissues; k++) {
	    	temp = 0.5 * ( hm + hp ) * (G[k].Evaluate(x+hp)+G[k].Evaluate(x-hm));
	        // MRF matrix fits number of tissues?
	        if( bMRF )
	        {
	        	temp = temp * MRFenergies[k] / denominatorMRF;
	        }
	        else
	        {
	        	temp = temp * _atlas.GetValue(k);
	        }

	        numerator[k] = temp;
	        denominator += temp;
	      }
	      delete[] MRFenergies;

	      for (k = 0; k < _number_of_tissues; k++) {
	        if (denominator != 0) {
	          double value = numerator[k]/denominator;
	          _output.SetValue(k, value);
	          if ((value < 0) || (value > 1)) {
	            cerr << "Probability value out of range = " << value << endl;
	            cerr << value << " " << k << " " << i << " " << _atlas.GetValue(k) <<  " " << _sigma[k] << endl;
	            exit(1);
	          }
	        } else {
	          cerr << "Division by 0 while computing probabilities" << endl;
	          cerr<<"tissue="<<k;
	          if (k==0) _output.SetValue(k, 1);
	          else     _output.SetValue(k, 0);
	        }
	      }
	    } else {
	      for (k = 0; k < _number_of_tissues - 1; k++) {
	        _output.SetValue(k, 0);
	      }
	      _output.SetValue(_number_of_tissues - 1, 1);
	    }
	    ptr++;
	    pm++;
	    _atlas.Next();
	    _output.Next();
	  }
	  delete[] numerator;
	  delete[] G;
}


void irtkEMClassification2ndOrderMRF::EStepMRF()
{
	  cout << "E-Step with 1st order MRF" << endl;
	  int i, k;
	  double x;
	  irtkGaussian* G = new irtkGaussian[_number_of_tissues];

	  for (k = 0; k < _number_of_tissues; k++) {
	    G[k].Initialise( _mi[k], _sigma[k]);
	  }

	  _atlas.First();
	  _output.First();
	  irtkRealPixel *ptr = _input.GetPointerToVoxels();
	  irtkRealPixel *pm = _mask.GetPointerToVoxels();

	  int per = 0;
	  double* numerator = new double[_number_of_tissues];
	  double denominator=0, temp=0;
	  bool bMRF = _number_of_tissues == _connectivity.Rows();
	  if(!bMRF)
	  {
		  cout << "Warning: number of tissues does not match size of connectivity matrix!" << endl;
	  }

	  for (i=0; i< _input.GetNumberOfVoxels(); i++) {
		if (i*10.0/_input.GetNumberOfVoxels() > per) {
		  per++;
		  cerr<<per<<"0%...";
		}
		denominator = 0;
		temp = 0;

	    if (*pm == 1) {
	      x = *ptr;
		  double hp = 0.5;
		  double hm = 0.5;
		  if( _isLogTransformed )
		  {
				hp = log( exp(x)+0.5 ) - x;
				hm = x - log( (exp(x) - 0.5) > 0 ? exp(x) - 0.5 : exp(x) );
		  }
	      double* MRFenergies = new double[_number_of_tissues];
	      double denominatorMRF = .0;

	      for (k = 0; k < _number_of_tissues; k++) {
	    	  MRFenergies[k] = _atlas.GetValue(k) * getMRFenergy(i,k);
	    	  denominatorMRF += MRFenergies[k];
	      }

	      for (k = 0; k < _number_of_tissues; k++) {
	    	temp = 0.5 * ( hm + hp ) * (G[k].Evaluate(x+hp)+G[k].Evaluate(x-hm));
	        //temp = G[k].Evaluate(x);
	        // MRF matrix fits number of tissues?
	        if( bMRF )
	        {
	        	temp = temp * MRFenergies[k] / denominatorMRF;
	        }
	        else
	        {
	        	temp = temp * _atlas.GetValue(k);
	        }

	        numerator[k] = temp;
	        denominator += temp;
	      }
	      delete[] MRFenergies;

	      for (k = 0; k < _number_of_tissues; k++) {
	        if (denominator != 0) {
	          double value = numerator[k]/denominator;
	          _output.SetValue(k, value);
	          if ((value < 0) || (value > 1)) {
	            cerr << "Probability value out of range = " << value << endl;
	            cerr << value << " " << k << " " << i << " " << _atlas.GetValue(k) <<  " " << _sigma[k] << endl;
	            exit(1);
	          }
	        } else {
	          cerr << "Division by 0 while computing probabilities" << endl;
	          cerr<<"tissue="<<k;
	          if (k==0) _output.SetValue(k, 1);
	          else     _output.SetValue(k, 0);
	        }
	      }
	    } else {
	      for (k = 0; k < _number_of_tissues - 1; k++) {
	        _output.SetValue(k, 0);
	      }
	      _output.SetValue(_number_of_tissues - 1, 1);
	    }
	    pm++;
	    ptr++;
	    _atlas.Next();
	    _output.Next();
	  }
	  delete[] numerator;
	  delete[] G;
}

void irtkEMClassification2ndOrderMRF::SetBiasField(irtkBiasField *biasfield)
{
	_biasfield = biasfield;
}

void irtkEMClassification2ndOrderMRF::GetBiasCorrectedImage(irtkRealImage &image)
{
  image = _input;
  //_biascorrection.Apply(_input);
}

void irtkEMClassification2ndOrderMRF::GetBiasField(irtkRealImage &image)
{
	irtkRealPixel *pm = _mask.GetPointerToVoxels();
	irtkRealPixel *ptrA = _input.GetPointerToVoxels();
	irtkRealPixel *ptrB = _uncorrected.GetPointerToVoxels();
	irtkRealPixel *ptrOutput = image.GetPointerToVoxels();

	for( int i = 0; i < image.GetNumberOfVoxels(); ++i )
	{
		if( *pm == 1 )
		{
			*ptrOutput = *ptrB - *ptrA;
		}
		else
		{
			*ptrOutput = 0;
		}
		ptrOutput++;
		ptrA++;
		ptrB++;
		pm++;
	}
}

bool irtkEMClassification2ndOrderMRF::isPVclass(int pvclass)
{
	if( pv_classes.find(pvclass) != pv_classes.end() ) return true;
	return false;
}

// Alternative cardoso implementation, not necessarily better!!
void irtkEMClassification2ndOrderMRF::MStepPV()
{
  int i, k;
  vector<double> mi_num(_number_of_tissues);
  vector<double> sigma_num(_number_of_tissues);
  vector<double> denom(_number_of_tissues);
  vector<bool> isPV(_number_of_tissues);

  for( int k = 0; k < _number_of_tissues; ++k )
  {
	isPV[k] = isPVclass(k);
  }

  for (k = 0; k < _number_of_tissues; k++) {
    mi_num[k] = 0;
  }

  for (k = 0; k < _number_of_tissues; k++) {
    sigma_num[k] = 0;
  }

  for (k = 0; k < _number_of_tissues; k++) {
    denom[k] = 0;
  }

  _output.First();
  irtkRealPixel *ptr = _input.GetPointerToVoxels();
  irtkRealPixel *pm = _mask.GetPointerToVoxels();

  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
    // Check for backgound
    if (*pm == 1) {
      for (k = 0; k < _number_of_tissues; k++) {
    	  if( !isPV[k] )
    	  {
    		  mi_num[k] += _output.GetValue(k) * *ptr;
    		  denom[k]  += _output.GetValue(k);
    	  }
      }
    }
    ptr++;
    pm++;
    _output.Next();
  }

  for (k = 0; k < _number_of_tissues; k++) {
	  if( !isPV[k] )
	  {
		if (denom[k] != 0) {
		  _mi[k] = mi_num[k] / denom[k];
		} else {
		  cerr << "Division by zero while computing tissue mean!" << endl;
		  exit(1);
		}
	  }
  }

  _output.First();
  ptr = _input.GetPointerToVoxels();
  pm = _mask.GetPointerToVoxels();

  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
    // Check for backgound
    if (*pm == 1) {
      for (k = 0; k <_number_of_tissues; k++) {
    	  if( !isPV[k] )
    	  {
    		  sigma_num[k] += (_output.GetValue(k) * (*ptr - _mi[k]) * (*ptr - _mi[k]));
    	  }
      }
    }
    ptr++;
    pm++;
    _output.Next();
  }

  for (k = 0; k <_number_of_tissues; k++) {
	  if( !isPV[k] )
	  {
		  cout << "sigma_num[" << k << "]=" << sigma_num[k] << endl;
		  cout << "denom[" << k << "]=" << denom[k] << endl;
		  _sigma[k] = sigma_num[k] / denom[k];
		  _sigma[k] = max( _sigma[k], 0.01 );
	  }
  }

  map<int,int>::iterator iter = pv_classes.begin();
  while( iter != pv_classes.end() )
  {
	  int position = iter->second;
	  int classA = pv_connections[position].first;
	  int classB = pv_connections[position].second;
	  int classPV = iter->first;
	  _mi[classPV] = (_mi[classA]+_mi[classB])/(float)(2.0);
	  _sigma[classPV] = sqrt(0.5*0.5*pow(_sigma[classA],2)+(0.5*0.5)*powf(_sigma[classB],2));
	  if(_sigma[classPV]<0.003) _sigma[classPV] = 0.003;
	  iter++;
  }
}

double irtkEMClassification2ndOrderMRF::LogLikelihood()
{
  int i, k;
  double temp, f;
  double x;
  cerr<< "Log likelihood: ";
  irtkGaussian* G = new irtkGaussian[_number_of_tissues];
  double* gv = new double[_number_of_tissues];

  for (k = 0; k < _number_of_tissues; k++) {
    G[k].Initialise( _mi[k], _sigma[k]);
  }

  irtkRealPixel *ptr = _input.GetPointerToVoxels();
  irtkRealPixel *pm = _mask.GetPointerToVoxels();
  _output.First();
  f = 0;
  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
    if (*pm == 1) {
      temp = 0;
      x = *ptr;

      for (k=0; k < _number_of_tissues; k++) {
		double hp = 0.5;
		double hm = 0.5;
		if( _isLogTransformed )
		{
			hp = log( exp(x)+0.5 ) - x;
			hm = x - log( (exp(x) - 0.5) > 0 ? exp(x) - 0.5 : exp(x) );
		}
        // Estimation of gaussian probability of intensity (*ptr) for tissue k
        gv[k] = 0.5 * ( hm + hp ) * (G[k].Evaluate(x+hp)+G[k].Evaluate(x-hm));

        // Probability that current voxel is from tissue k
        temp += gv[k] * _output.GetValue(k);

      }
      if( temp < 0 || temp > 1)
      {
          cerr << "Could not compute likelihood, probability out of range = " << temp << endl;
          exit(1);
      }

      f += log(temp);
    }
    ptr++;
    pm++;
    _output.Next();
  }

  f = -f;
  double diff, rel_diff;
  diff = _f-f;

  if (_f == 0) rel_diff = 1;
  else rel_diff = diff/_f;

  _f=f;

  cerr << "f= "<< f << " diff = " << diff << " rel_diff = " << rel_diff <<endl;
  delete[] G;
  delete[] gv;

  return rel_diff;
}

void irtkEMClassification2ndOrderMRF::EStep()
{
  int i, k;
  double x;

  irtkGaussian* G = new irtkGaussian[_number_of_tissues];

  for (k = 0; k < _number_of_tissues; k++) {
    G[k].Initialise( _mi[k], _sigma[k]);
  }

  _atlas.First();
  _output.First();
  irtkRealPixel *ptr = _input.GetPointerToVoxels();
  irtkRealPixel *pm = _mask.GetPointerToVoxels();
  for (i=0; i< _input.GetNumberOfVoxels(); i++) {
    double* gv = new double[_number_of_tissues];
    double* numerator = new double[_number_of_tissues];
    double denominator=0, temp=0;
    if (*pm == 1) {
      x = *ptr;
	  double hp = 0.5;
	  double hm = 0.5;
	  if( _isLogTransformed )
	  {
			hp = log( exp(x)+0.5 ) - x;
			hm = x - log( (exp(x) - 0.5) > 0 ? exp(x) - 0.5 : exp(x) );
	  }
      for (k = 0; k < _number_of_tissues; k++) {
        temp = 0.5 * ( hm + hp ) * (G[k].Evaluate(x+hp)+G[k].Evaluate(x-hm));
    	//temp = G[k].Evaluate(x);
        gv[k] = temp;
        temp = temp * _atlas.GetValue(k);
        numerator[k] = temp;
        denominator += temp;
      }
      for (k = 0; k < _number_of_tissues; k++) {
        if (denominator != 0) {
          double value = numerator[k]/denominator;
          _output.SetValue(k, value);
          if ((value < 0) || (value > 1)) {
            cerr << "Probability value out of range = " << value << endl;
            cerr << value << " " << k << " " << i << " " << _atlas.GetValue(k) <<  " " << _sigma[k] << endl;
            exit(1);
          }
        } else {
          cerr << "Division by 0 while computing probabilities" << endl;
          cerr<<"tissue="<<k;
          if (k==0) _output.SetValue(k, 1);
          else     _output.SetValue(k, 0);
        }
      }
    } else {
      for (k = 0; k < _number_of_tissues - 1; k++) {
        _output.SetValue(k, 0);
      }
      _output.SetValue(_number_of_tissues - 1, 1);
    }
    ptr++;
    pm++;
    _atlas.Next();
    _output.Next();
    delete[] gv;
    delete[] numerator;
  }
  delete[] G;

}

