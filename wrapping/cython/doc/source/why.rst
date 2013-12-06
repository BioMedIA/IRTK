Why Python?
===========

.. code-block:: python

    ls *.nii | python -c 'import sys, irtk; [sys.stdout.write( str(irtk.imread(line.rstrip(), dtype="float32").max())+"\n") for line in sys.stdin]'

    
Why another Python libraries for medical imaging?
=================================================

.. code-block:: python
                
    import itk

    pixelType = itk.UC
    imageType = itk.Image[pixelType, 3]
    readerType = itk.ImageFileReader[imageType]
    writerType = itk.ImageFileWriter[imageType]
    reader = readerType.New()
    reader.SetFileName( "input.nii" )
    reader.Update()

    itk2np = itk.PyBuffer[imageType]
    data = itk2np.GetArrayFromImage( reader.GetOutput() )

    ...


.. code-block:: python
                
    import SimpleITK as sitk

    img = sitk.ReadImage( "input.nii" )
    data = sitk.GetArrayFromImage( img )

    ...

    output = sitk.GetImageFromArray( data )
    output.SetSpacing( img.GetSpacing() )
    output.SetOrigin( img.GetOrigin() )
    output.SetDirection( img.GetDirection() )
    sitk.WriteImage( output, "output.nii" ) 

.. code-block:: python
                
    import irtk

    img = irtk.imread( "input.nii" )

    ...

    irtk.imwrite( "output.nii", img )

