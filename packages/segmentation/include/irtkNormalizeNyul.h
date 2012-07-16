#ifndef IRTKNORMALIZENYUL_H_
#define IRTKNORMALIZENYUL_H_

#include <irtkImage.h>
#include <irtkRegistration.h>

class irtkNormalizeNyul{

private:
	irtkRealImage _target;
	irtkRealImage _source;
	int _source_padding;
	int _target_padding;



public:
	irtkNormalizeNyul(irtkRealImage source, irtkRealImage target);
	void SetMask(irtkRealImage source_mask, irtkRealImage target_mask);
	void SetPadding(int source_padding, int target_padding);
	void Run();
	irtkRealImage GetOutput();
	irtkRealImage GetTarget();
};

#endif
