/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include "irtkSegmentationFunction.h"
#ifdef HAS_OPENCV
/* these settings affect the quality of detection: change with care */
void on_mouse( int event, int x, int y, int flags, void* param )
{
    switch( event )
    {
    case CV_EVENT_LBUTTONDOWN:
        {
            ((CvRect*)param)->x = x; ((CvRect*)param)->y = y;
        }
        break;
    case CV_EVENT_LBUTTONUP:
        {
            ((CvRect*)param)->width = x - ((CvRect*)param)->x;
            ((CvRect*)param)->height = y - ((CvRect*)param)->y;
        }
        break;
    }
}
/// detect possible presence of the object at threshold == 3
irtkAffineTransformation irtkSegmentationFunction::DetectObject (irtkRealImage *threshold, irtkGreyImage *target, CvHaarClassifierCascade *classifier, int oxsize, double fscale, int size)
{

  int i,k;
  //irtkPoint* heart_center = new irtkPoint[4];
  irtkAffineTransformation transformation;
  if(classifier == NULL) {
    cout << "No Classifier given"<<endl;
    exit(1);
  } else {
    cout << "Evaluating Object Center"<<endl;
    irtkImageToOpenCv<irtkGreyPixel> itocg;
    irtkImageToOpenCv<irtkRealPixel> itocr;
    //target->GetMinMax(&min,&max);
    double scale = oxsize/target->GetXSize();
    irtkAffineTransformation tmptransformation[20];
    //create mem storage for detection
    static CvMemStorage* storage = 0;
    storage = cvCreateMemStorage(0);
    cvClearMemStorage( storage );

    //create target image
    IplImage* pimage= cvCreateImage(cvSize(target->GetX(), target->GetY()), IPL_DEPTH_8U, 1);
    //create threshold image
    IplImage* pths= cvCreateImage(cvSize(threshold->GetX(), threshold->GetY()), IPL_DEPTH_8U, 1);
    //create small_img and convert to gray images
    IplImage* small_img = cvCreateImage( cvSize( cvRound (pimage->width/scale),
                                         cvRound (pimage->height/scale)),
                                         IPL_DEPTH_8U, 1 );
    IplImage* small_ths = cvCreateImage( cvSize( cvRound (pths->width/scale),
                                         cvRound (pths->height/scale)),
                                         IPL_DEPTH_8U, 1 );
    CvPoint center;
    CvRect* r;
    CvPoint textcenter;
    int count = 0;
    // find center of z axis
    for (k = (target->GetZ() - 1)/ 2 - round(target->GetZ()/8.0); k < (target->GetZ() - 1)/ 2 +round(target->GetZ()/8.0+1); k++) {
      //write pixel
      itocg.SetInput(target);
      itocg.SetOutput(pimage);
      itocg.Run(k);
      /*for (j = 0; j < target->GetY(); j++) {
      for (i = 0; i < target->GetX(); i++) {
      //input[i](x, y, z) = inputsequence(x,y,z,i);
      int tmp = (target->GetAsDouble(i,j,k,0) * 256 /max);
      pimage->imageData[j*pimage->widthStep + i] = tmp;
    }
    }*/
      //threshold->Write("threshold.gipl");
      //write pixel
      itocr.SetInput(threshold);
      itocr.SetOutput(pths);
      itocr.Run(k);
      /*for (j = 0; j < threshold->GetY(); j++) {
      for (i = 0; i < threshold->GetX(); i++) {
      //input[i](x, y, z) = inputsequence(x,y,z,i);
      int tmp = threshold->GetAsDouble(i,j,k,0);
      pths->imageData[j*pths->widthStep + i] = tmp;
    }
    }*/

      //cvCvtColor( pimage, gray, CV_BGR2GRAY );
      //cvResize( gray, small_img, CV_INTER_LINEAR );
      cvResize( pimage, small_img, CV_INTER_LINEAR );
      //cvSaveImage("pths.jpg",pths);
      cvResize( pths, small_ths, CV_INTER_LINEAR );
      cvEqualizeHist( small_img, small_img );
      //cvSaveImage("small_img.jpg",small_img);
      //cvSaveImage("small_ths.jpg",small_ths);

      //detect interest region using modified haar detect object function
      CvSeq* heart = cvTHaarDetectObjects( small_img, small_ths, classifier, storage,
                                           1.05, 2, 0, cvSize(size,size) );//1 not canny 0 canny
      //CvSeq* heart = cvHaarDetectObjects( small_img, classifier, storage,
      //               1.1, 2, 1, cvSize(50,50) );

      cout << "Evaluation Slice "<< k <<" Done"<<endl;

      //write interest region to affine transformation
      i = (heart ? heart->total : 0);
      if ( i > 1) {
        cout <<"multipule objects detected, we know this is wrong."<<endl;
        for( i = 0; i < heart->total; i++ ) {
          count++;
          r = (CvRect*)cvGetSeqElem( heart, i );
          center.x = cvRound((r->x + r->width*0.5)*scale);
          center.y = cvRound((r->y + r->height*0.5)*scale);
          //textcenter.x = center.x - 10;
          //textcenter.y = center.y + 10;
          //int  radius = cvRound((r->width + r->height)*0.25*scale);
          //cvCircle( pimage, center, radius, CV_RGB(200,200,200), 3, 8, 0 );
          //sprintf(numberbuff, "%d",i);
          //cvPutText( pimage, numberbuff, textcenter, &font, CV_RGB(200,200,200));
          //this is for test only
          cout << "Object Center Location: " << center.x << " " << center.y << endl;
          cout << "Object Box Scale: " << r->width*scale << " " << r->height*scale << endl;
          tmptransformation[count-1].Put(0,center.x);
          tmptransformation[count-1].Put(1,center.y);
          tmptransformation[count-1].Put(6,r->width*target->GetXSize()*scale);
          tmptransformation[count-1].Put(7,r->height*target->GetYSize()*scale);
          if (count >= 20)
            break;
          //the above is for test only
        }
        /*
        i = heart->total;
        while(i < 0 || i >= heart->total){
        cvNamedWindow( "result", 1 );
        cvShowImage( "result", pimage );
        cout<<"Please press enter to close the window \nand select the correct one by entering the number: ";
        cvWaitKey();
        cvDestroyWindow("result");
        cin>>i;
      }
        count++;
        */
        //delete []numberbuff;
      } else if(i==1) {
        i = 0;
        count++;
        //this is for test only
        r = (CvRect*)cvGetSeqElem( heart, i );
        center.x = cvRound((r->x + r->width*0.5)*scale);
        center.y = cvRound((r->y + r->height*0.5)*scale);
        cout << "Object Center Location: " << center.x << " " << center.y << endl;
        cout << "Object Box Scale: " << r->width*scale << " " << r->height*scale << endl;
        tmptransformation[count-1].Put(0,center.x);
        tmptransformation[count-1].Put(1,center.y);
        tmptransformation[count-1].Put(6,r->width*target->GetXSize()*scale);
        tmptransformation[count-1].Put(7,r->height*target->GetYSize()*scale);
        if (count >= 20)
          break;
        //the above is for test only
      } else {
        cerr <<"no object detected, this is wrong!"<<endl;
        continue;
      }
      if (count >= 20)
        break;
      /*r = (CvRect*)cvGetSeqElem( heart, i );
      center.x = cvRound((r->x + r->width*0.5)*scale);
      center.y = cvRound((r->y + r->height*0.5)*scale);
      cout << "Evaluation Done"<<endl;
      cout << "Heart Center Location: "<< count << " " << center.x << " " << center.y << endl;
      cout << "Heart Box Scale: " << r->width << " " << r->height << endl;
      tmptransformation[count-1].Put(0,center.x);
      tmptransformation[count-1].Put(1,center.y);
      tmptransformation[count-1].Put(6,r->width*target->GetXSize());
      tmptransformation[count-1].Put(7,r->height*target->GetYSize());*/

    }
    if (count == 0) {
      cerr <<"no object detected at all, this is really wrong!"<<endl;
      cout <<"Sorry but Please select the interest region by draging the mouse"<<endl;
      r = (CvRect*)cvAlloc(sizeof(CvRect));
      CvFont font;
      cvInitFont( &font, CV_FONT_VECTOR0,0.3, 0.3);
      char numberbuff[255];
      //create pimage
      itocg.Run(k);
      /*for (j = 0; j < target->GetY(); j++) {
      for (i = 0; i < target->GetX(); i++) {
      int tmp = (target->GetAsDouble(i,j,k-1,0) * 256 /max);
      pimage->imageData[j*pimage->widthStep + i] = tmp;
    }
    }*/

      cvNamedWindow( "result", 1 );
      //create event
      cvSetMouseCallback( "result", on_mouse, r);
      while(!count) {
        itocg.Run(k);
        /*
        for (j = 0; j < target->GetY(); j++) {
        for (i = 0; i < target->GetX(); i++) {
        int tmp = (target->GetAsDouble(i,j,k-1,0) * 256 /max);
        pimage->imageData[j*pimage->widthStep + i] = tmp;
      }
      }*/
        sprintf(numberbuff, "Drag mouse to select region");
        textcenter.x = 10;
        textcenter.y = pimage->height - 40;
        cvPutText( pimage, numberbuff, textcenter, &font, CV_RGB(200,200,200));

        sprintf(numberbuff, "Press enter to continue");
        textcenter.x = 10;
        textcenter.y = pimage->height - 20;
        cvPutText( pimage, numberbuff, textcenter, &font, CV_RGB(200,200,200));

        cvShowImage( "result", pimage );
        //hint
        cvWaitKey();
        itocg.Run(k);
        /*for (j = 0; j < target->GetY(); j++) {
        for (i = 0; i < target->GetX(); i++) {
        int tmp = (target->GetAsDouble(i,j,k-1,0) * 256 /max);
        pimage->imageData[j*pimage->widthStep + i] = tmp;
      }
      }*/
        //final check
        center.x = cvRound((r->x + r->width*0.5));
        center.y = cvRound((r->y + r->height*0.5));
        CvPoint p1= cvPoint(r->x,r->y);
        CvPoint p2= cvPoint(r->x+r->width,r->y+r->height);
        cvRectangle( pimage, p1, p2, CV_RGB(200,200,200), 3, 8, 0 );
        sprintf(numberbuff, "Enter 1 to confirm else to repick");
        textcenter.x = 10;
        textcenter.y = r->height + r->y + 10;
        cvPutText( pimage, numberbuff, textcenter, &font, CV_RGB(200,200,200));
        cvShowImage( "result", pimage );
        count = ((char)cvWaitKey() == '1');
      }
      cout << "Object Center Location: " << center.x << " " << center.y << endl;
      cout << "Object Box Scale: " << r->width << " " << r->height << endl;
      transformation.Put(0,center.x);
      transformation.Put(1,center.y);
      //let the bounding box larger
      transformation.Put(6,r->width*target->GetXSize()*1.1);
      transformation.Put(7,r->height*target->GetYSize()*1.1);
      cvFree(&r);
    } else if(count == 1) {
      cerr<<"only one object detected in three slice, doubt the accuracy!"<<endl;
      transformation.Put(0,tmptransformation[0].Get(0));
      transformation.Put(1,tmptransformation[0].Get(1));
      //let the bounding box larger
      transformation.Put(6,tmptransformation[0].Get(6)*1.1);
      transformation.Put(7,tmptransformation[0].Get(7)*1.1);
    } else {
      cout<<"multiple objects detected combining..."<<endl;
      //find outlayer and kill them first

      //now combine
      while(count > 2) {
        //evaluate average
        double mx = 0,my = 0;
        int exclude = 0;
        double dmax = 0, tmax;
        for( k = 0 ; k < count; k++) {
          mx += tmptransformation[k].Get(0);
          my += tmptransformation[k].Get(1);
        }
        mx = mx/count;
        my = my/count;
        //exclude the most distant one
        for( k = 0 ; k < count; k++) {
          tmax = (tmptransformation[k].Get(0) - mx)*(tmptransformation[k].Get(0) - mx)
                 + (tmptransformation[k].Get(1) - my)*(tmptransformation[k].Get(1) - my);
          if(tmax>dmax) {
            dmax = tmax;
            exclude = k;
          }
        }
        tmptransformation[exclude].Put(0,tmptransformation[count - 1].Get(0));
        tmptransformation[exclude].Put(1,tmptransformation[count - 1].Get(1));
        tmptransformation[exclude].Put(6,tmptransformation[count - 1].Get(6));
        tmptransformation[exclude].Put(7,tmptransformation[count - 1].Get(7));
        count--;
      }
      //combine the rest two
      transformation.Put(0,(tmptransformation[0].Get(0)+tmptransformation[1].Get(0))/2);
      transformation.Put(1,(tmptransformation[0].Get(1)+tmptransformation[1].Get(1))/2);
      //let the bounding box larger
      transformation.Put(6,(tmptransformation[0].Get(6)>tmptransformation[1].Get(6)
                            ?tmptransformation[0].Get(6):tmptransformation[1].Get(6)) * fscale);
      transformation.Put(7,(tmptransformation[0].Get(7)>tmptransformation[1].Get(7)
                            ?tmptransformation[0].Get(7):tmptransformation[1].Get(7)) * fscale);
      cout << "Final Evaluation Done"<<endl;
      cout << "Object Center Location generated from "<< count << " instances:" << endl;
      cout << "Object Box Center & Scale: " << transformation.Get(0) << " " << transformation.Get(1) <<" "<< transformation.Get(6) << "mm " << transformation.Get(7) << "mm"<< endl;
    }
    // release images
    cvReleaseImage(&pimage);
    cvReleaseImage(&pths);
    cvReleaseImage(&small_img);
    cvReleaseImage(&small_ths);
    cvReleaseMemStorage(&storage);
  }
  return transformation;
}

irtkAffineTransformation irtkSegmentationFunction::DetectObject (irtkGreyImage *target, CvHaarClassifierCascade *classifier, int oxsize, double fscale, int size)
{

  int i,k;
  //irtkPoint* heart_center = new irtkPoint[4];
  irtkAffineTransformation transformation;
  if(classifier == NULL) {
    cerr << "do not have classifier" <<endl;
    exit(1);
  }
  cout << "Evaluating Object Center"<<endl;
  irtkImageToOpenCv<irtkGreyPixel> itocg;
  irtkImageToOpenCv<irtkRealPixel> itocr;
  //target->GetMinMax(&min,&max);
  double scale = oxsize/target->GetXSize();
  irtkAffineTransformation tmptransformation[20];
  //create mem storage for detection
  static CvMemStorage* storage = 0;
  storage = cvCreateMemStorage(0);
  cvClearMemStorage( storage );

  //create target image
  IplImage* pimage= cvCreateImage(cvSize(target->GetX(), target->GetY()), IPL_DEPTH_8U, 1);
  //create small_img and convert to gray images
  IplImage* small_img = cvCreateImage( cvSize( cvRound (pimage->width/scale),
                                       cvRound (pimage->height/scale)),
                                       IPL_DEPTH_8U, 1 );
  CvPoint center;
  CvRect* r;
  CvPoint textcenter;
  int count = 0;
  // find center of z axis
  for (k = (target->GetZ() - 1)/ 2 - round(target->GetZ()/8.0); k < (target->GetZ() - 1)/ 2 +round(target->GetZ()/8.0+1); k++) {
    //write pixel
    itocg.SetInput(target);
    itocg.SetOutput(pimage);
    itocg.Run(k);
    cvResize( pimage, small_img, CV_INTER_LINEAR );
    cvEqualizeHist( small_img, small_img );

    //detect interest region using modified haar detect object function
    CvSeq* heart = cvHaarDetectObjects( small_img, classifier, storage,
                                        1.05, 2, 0, cvSize(size,size) );
    cout << "Evaluation Slice "<< k <<" Done"<<endl;

    //write interest region to affine transformation
    i = (heart ? heart->total : 0);
    if ( i > 1) {
      cout <<"multipule objects detected, we know this is wrong."<<endl;
      for( i = 0; i < heart->total; i++ ) {
        count++;
        r = (CvRect*)cvGetSeqElem( heart, i );
        center.x = cvRound((r->x + r->width*0.5)*scale);
        center.y = cvRound((r->y + r->height*0.5)*scale);
        cout << "Object Center Location: " << center.x << " " << center.y << endl;
        cout << "Object Box Scale: " << r->width*scale << " " << r->height*scale << endl;
        tmptransformation[count-1].Put(0,center.x);
        tmptransformation[count-1].Put(1,center.y);
        tmptransformation[count-1].Put(6,r->width*target->GetXSize()*scale);
        tmptransformation[count-1].Put(7,r->height*target->GetYSize()*scale);
        if (count >= 20)
          break;
        //the above is for test only
      }
    } else if(i==1) {
      i = 0;
      count++;
      //this is for test only
      r = (CvRect*)cvGetSeqElem( heart, i );
      center.x = cvRound((r->x + r->width*0.5)*scale);
      center.y = cvRound((r->y + r->height*0.5)*scale);
      cout << "Object Center Location: " << center.x << " " << center.y << endl;
      cout << "Object Box Scale: " << r->width*scale << " " << r->height*scale << endl;
      tmptransformation[count-1].Put(0,center.x);
      tmptransformation[count-1].Put(1,center.y);
      tmptransformation[count-1].Put(6,r->width*target->GetXSize()*scale);
      tmptransformation[count-1].Put(7,r->height*target->GetYSize()*scale);
      if (count >= 20)
        break;
      //the above is for test only
    } else {
      cerr <<"no object detected, this is wrong!"<<endl;
      continue;
    }
    if (count >= 20)
      break;

  }
  if (count == 0) {
    cerr <<"no object detected at all, this is really wrong!"<<endl;
    cout <<"Sorry but Please select the interest region by draging the mouse"<<endl;
    r = (CvRect*)cvAlloc(sizeof(CvRect));
    CvFont font;
    cvInitFont( &font, CV_FONT_VECTOR0,0.3, 0.3);
    char numberbuff[255];
    //create pimage
    itocg.Run(k);
    cvNamedWindow( "result", 1 );
    //create event
    cvSetMouseCallback( "result", on_mouse, r);
    while(!count) {
      itocg.Run(k);
      sprintf(numberbuff, "Drag mouse to select region");
      textcenter.x = 10;
      textcenter.y = pimage->height - 40;
      cvPutText( pimage, numberbuff, textcenter, &font, CV_RGB(200,200,200));

      sprintf(numberbuff, "Press enter to continue");
      textcenter.x = 10;
      textcenter.y = pimage->height - 20;
      cvPutText( pimage, numberbuff, textcenter, &font, CV_RGB(200,200,200));

      cvShowImage( "result", pimage );
      //hint
      cvWaitKey();
      itocg.Run(k);
      //final check
      center.x = cvRound((r->x + r->width*0.5));
      center.y = cvRound((r->y + r->height*0.5));
      CvPoint p1= cvPoint(r->x,r->y);
      CvPoint p2= cvPoint(r->x+r->width,r->y+r->height);
      cvRectangle( pimage, p1, p2, CV_RGB(200,200,200), 3, 8, 0 );
      sprintf(numberbuff, "Enter 1 to confirm else to repick");
      textcenter.x = 10;
      textcenter.y = r->height + r->y + 10;
      cvPutText( pimage, numberbuff, textcenter, &font, CV_RGB(200,200,200));
      cvShowImage( "result", pimage );
      count = ((char)cvWaitKey() == '1');
    }
    cout << "Object Center Location: " << center.x << " " << center.y << endl;
    cout << "Object Box Scale: " << r->width << " " << r->height << endl;
    transformation.Put(0,center.x);
    transformation.Put(1,center.y);
    //let the bounding box larger
    transformation.Put(6,r->width*target->GetXSize()*1.1);
    transformation.Put(7,r->height*target->GetYSize()*1.1);
    cvFree(&r);
  } else if(count == 1) {
    cerr<<"only one object detected in three slice, doubt the accuracy!"<<endl;
    transformation.Put(0,tmptransformation[0].Get(0));
    transformation.Put(1,tmptransformation[0].Get(1));
    //let the bounding box larger
    transformation.Put(6,tmptransformation[0].Get(6)*1.1);
    transformation.Put(7,tmptransformation[0].Get(7)*1.1);
  } else {
    cout<<"multiple objects detected combining..."<<endl;
    //find outlayer and kill them first

    //now combine
    while(count > 2) {
      //evaluate average
      double mx = 0,my = 0;
      int exclude = 0;
      double dmax = 0, tmax;
      for( k = 0 ; k < count; k++) {
        mx += tmptransformation[k].Get(0);
        my += tmptransformation[k].Get(1);
      }
      mx = mx/count;
      my = my/count;
      //exclude the most distant one
      for( k = 0 ; k < count; k++) {
        tmax = (tmptransformation[k].Get(0) - mx)*(tmptransformation[k].Get(0) - mx)
               + (tmptransformation[k].Get(1) - my)*(tmptransformation[k].Get(1) - my);
        if(tmax>dmax) {
          dmax = tmax;
          exclude = k;
        }
      }
      tmptransformation[exclude].Put(0,tmptransformation[count - 1].Get(0));
      tmptransformation[exclude].Put(1,tmptransformation[count - 1].Get(1));
      tmptransformation[exclude].Put(6,tmptransformation[count - 1].Get(6));
      tmptransformation[exclude].Put(7,tmptransformation[count - 1].Get(7));
      count--;
    }
    //combine the rest two
    transformation.Put(0,(tmptransformation[0].Get(0)+tmptransformation[1].Get(0))/2);
    transformation.Put(1,(tmptransformation[0].Get(1)+tmptransformation[1].Get(1))/2);
    //let the bounding box larger
    transformation.Put(6,(tmptransformation[0].Get(6)>tmptransformation[1].Get(6)
                          ?tmptransformation[0].Get(6):tmptransformation[1].Get(6)) * fscale);
    transformation.Put(7,(tmptransformation[0].Get(7)>tmptransformation[1].Get(7)
                          ?tmptransformation[0].Get(7):tmptransformation[1].Get(7)) * fscale);
    cout << "Final Evaluation Done"<<endl;
    cout << "Object Center Location generated from "<< count << " instances:" << endl;
    cout << "Object Box Center & Scale: " << transformation.Get(0) << " " << transformation.Get(1) <<" "<< transformation.Get(6) << "mm " << transformation.Get(7) << "mm"<< endl;
  }
  // release images
  cvReleaseImage(&pimage);
  cvReleaseImage(&small_img);
  cvReleaseMemStorage(&storage);
  return transformation;
}

#else
irtkAffineTransformation irtkSegmentationFunction::DetectObject (irtkRealImage *, irtkGreyImage *, int, double)
{
  cerr << "irtkSegmentationFunction::DetectObject: Must have OpenCV enabled in order to detect an object" << endl;
  exit(1);
}
#endif

void irtkSegmentationFunction::GenerateBox(irtkAffineTransformation& interestregion, irtkGreyImage* interest, irtkPoint& ip1, irtkPoint& ip2, irtkGreyImage* utarget)
{
  int i,j,k;
  for (k = 0; k < utarget->GetZ(); k++) {
    for (j = 0; j < utarget->GetY(); j++) {
      for (i = 0; i < utarget->GetX(); i++) {
        interest->PutAsDouble(i,j,k,0);
      }
    }
  }

  ip1._x = interestregion.Get(0) - interestregion.Get(6)/utarget->GetXSize()/2;
  ip1._y = interestregion.Get(1) - interestregion.Get(7)/utarget->GetYSize()/2;
  ip1._z = 0;

  ip2._x = interestregion.Get(0) + interestregion.Get(6)/utarget->GetXSize()/2;
  ip2._y = interestregion.Get(1) + interestregion.Get(7)/utarget->GetYSize()/2;
  ip2._z = interest->GetZ();

  if (ip1._x < 0) ip1._x = 0;	if (ip1._x > interest->GetX()) ip1._x = interest->GetX();
  if (ip1._y < 0) ip1._y = 0;	if (ip1._y > interest->GetY()) ip1._y = interest->GetY();
  if (ip2._x < 0) ip2._x = 0;	if (ip2._x > interest->GetX()) ip2._x = interest->GetX();
  if (ip2._y < 0) ip2._y = 0;	if (ip2._y > interest->GetY()) ip2._y = interest->GetY();

  for (k = round(ip1._z); k < round(ip2._z); k++) {
    for (j = round(ip1._y); j < round(ip2._y); j++) {
      for (i = round(ip1._x); i < round(ip2._x); i++) {
        interest->PutAsDouble(i,j,k,2);
      }
    }
  }
}

void irtkSegmentationFunction::FixThreshold(irtkRealImage* threshold, irtkRealImage* region)
{
  int i,j,k;
  for(k=0;k<threshold->GetZ();k++) {
    for(j=0;j<threshold->GetY();j++) {
      for(i=0;i<threshold->GetX();i++) {
        if(threshold->GetAsDouble(i,j,k) != 2 && region->GetAsDouble(i,j,k) == 3) {
          threshold->PutAsDouble(i,j,k,3);
        }
      }
    }
  }
}

//region grow in image coordinate
void irtkSegmentationFunction::RegionGrow(irtkPoint& center, irtkGreyImage* target, irtkRealImage* threshold, int minT, int value)
{
  irtkPointSet region;
  int bsize,asize,i,j;
  int average;
  bsize = 0; asize = 1;
  region.Add(center);

  threshold->Initialize(target->GetImageAttributes());
  threshold->PutAsDouble(center._x,center._y,center._z,value);

  while(bsize != asize) {
    bsize = region.Size();
    //get average
    average = 0;
    for(i = 0; i < region.Size(); i++) {
      average += target->GetAsDouble(region(i)._x,region(i)._y,region(i)._z);
    }
    average = average/region.Size();
    for(i = 0; i < region.Size(); i++) {
      //regiongrowing here
      irtkPoint current,candidate[6];
      current = region(i);
      candidate[0]._x = region(i)._x-1;  candidate[0]._y = region(i)._y;  candidate[0]._z = region(i)._z;
      candidate[1]._x = region(i)._x+1;  candidate[1]._y = region(i)._y;  candidate[1]._z = region(i)._z;
      candidate[2]._x = region(i)._x;  candidate[2]._y = region(i)._y-1;  candidate[2]._z = region(i)._z;
      candidate[3]._x = region(i)._x;  candidate[3]._y = region(i)._y+1;  candidate[3]._z = region(i)._z;
      candidate[4]._x = region(i)._x;  candidate[4]._y = region(i)._y;  candidate[4]._z = region(i)._z-1;
      candidate[5]._x = region(i)._x;  candidate[5]._y = region(i)._y;  candidate[5]._z = region(i)._z+1;

      for(j = 0; j < 6; j++) {
        //validate candidate[i]
        if( candidate[j]._x > 0 && candidate[j]._x < target->GetX()
            && candidate[j]._y > 0 && candidate[j]._y < target->GetY()
            && candidate[j]._z > 0 && candidate[j]._z < target->GetZ()
            &&target->GetAsDouble(candidate[j]._x,candidate[j]._y,candidate[j]._z) >= average - minT
            &&(target->GetAsDouble(candidate[j]._x,candidate[j]._y,candidate[j]._z) <= average + minT || value == 2)
            && !threshold->GetAsDouble(candidate[j]._x,candidate[j]._y,candidate[j]._z)) {
          threshold->PutAsDouble(candidate[j]._x,candidate[j]._y,candidate[j]._z,value);
          region.Add(candidate[j]);
        }
      }
    }
    asize = region.Size();
  }

}

void irtkSegmentationFunction::EvaluateGraphCut (irtkGreyImage *threshold, irtkGreyImage **input,irtkRealImage **atlas, int numberofinput, int numberofatlas, double timeweight,int cutmode,int connectmode, double dataweight,int padding, int *c)
{
  double treshold = 0.0001;
  int iterations = 20;
  int i,n;
  irtkEMClassificationMultiComp *classification;
  irtkRealImage background;
  irtkRealImage **test;
  test = new irtkRealImage*[numberofatlas+1];

  cout << "Initialize classification"<<endl;
  classification= new irtkEMClassificationMultiComp(numberofatlas, atlas, NULL, c);
  cout << "Initialize finished"<<endl;
  classification->SetPadding(padding);
  classification->SetInput(*input[0]);
  classification->Initialise();

  double rel_diff,previous;
  cout << "Segmenting tissues: Background, Blood, Cardiac wall"<<endl;
  i=0; rel_diff = 100;
  do {
      previous = rel_diff;
      cout << "Iteration = " << i+1 << " / " << iterations << endl;
      rel_diff = classification->Iterate(i);
      i++;
  } while ((rel_diff>treshold)&&(i<iterations)&&(rel_diff < previous));

  classification->ConstructSegmentation(background);

  //weights for graph cut
  for(n = 0; n < numberofatlas + 1; n ++) {
      test[n] = new irtkRealImage(input[0]->GetImageAttributes());
      classification->GetProbMap(n,*test[n]);
  }
  delete classification;

  //add weights for graph cut
  irtkImageGraphCut<irtkGreyPixel> graphcut;
  graphcut.SetInput(numberofinput,input,numberofatlas+1,test);
  graphcut.SetOutput(threshold);
  graphcut.Setdt(timeweight);

  graphcut.SetMode(cutmode);
  graphcut.Run(dataweight,connectmode);

  //finalize
   for(n = 0; n < numberofatlas + 1; n ++) {
    delete test[n];
   }
  delete []test;
}

void irtkSegmentationFunction::EvaluateGraphCut (irtkRealImage **threshold, irtkRealImage **input,irtkRealImage **atlas, int numberofinput, int numberofatlas, double timeweight,int cutmode,int connectmode, double dataweight, int padding, int *c)
{
  double treshold = 0.0001;
  int iterations = 20;
  int i,n,m;
  irtkEMClassificationMultiComp *classification;
  irtkRealImage **tmpatlas;
  irtkRealImage *test;
  double *datacost,count,*imageoffset;
  irtkLinearInterpolateImageFunction *interpolator = new irtkLinearInterpolateImageFunction;
  test = new irtkRealImage[numberofatlas+1];
  tmpatlas = new irtkRealImage*[numberofatlas];
  imageoffset = new double[numberofinput];
  count = 0;
  for(n=0;n<numberofinput;n++) {
      imageoffset[n] = count;
      count += input[n]->GetNumberOfVoxels();
  }
  datacost = new double[round(count)*(numberofatlas+1)];

  for(n=0;n<numberofinput;n++) {
    cout << "Initialize classification"<<endl;
    for(i=0;i<numberofatlas;i++) {
      tmpatlas[i] = new irtkRealImage(input[n]->GetImageAttributes());
      // Create image transformation
      irtkImageTransformation *imagetransformation =
        new irtkImageTransformation;
      irtkTransformation *transformation = new irtkRigidTransformation;

      imagetransformation->SetInput(atlas[i],transformation);
      imagetransformation->SetOutput(tmpatlas[i]);
      imagetransformation->PutInterpolator(interpolator);

      // Transform image
      imagetransformation->Run();
      delete transformation;
      delete imagetransformation;
    }
    classification= new irtkEMClassificationMultiComp(numberofatlas, tmpatlas, NULL, c);
    cout << "Initialize finished"<<endl;
    classification->SetPadding(padding);
    classification->SetInput(*input[n]);
    classification->Initialise();

    double rel_diff,previous;
    cout << "Segmenting tissues: Background, Blood, Cardiac wall"<<endl;
    i=0; rel_diff = 100;
    do {
      previous = rel_diff;
      cout << "Iteration = " << i+1 << " / " << iterations << endl;
      rel_diff = classification->Iterate(i);
      i++;
    } while ((rel_diff>treshold)&&(i<iterations)&&(rel_diff < previous));

    //weights for graph cut
    for(m = 0; m < numberofatlas + 1; m ++) {
      classification->GetProbMap(m,test[m]);
      irtkRealPixel *ptr = test[m].GetPointerToVoxels();
      for(i = 0; i < test[m].GetNumberOfVoxels(); i++){
          datacost[round(imageoffset[n]+i)*(numberofatlas + 1)+m]=*ptr;
          ptr++;
      }
    }
    delete classification;
    for(i=0;i<numberofatlas;i++) {
      delete tmpatlas[i];
    }
  }

  //add weights for graph cut
  irtkMultiImageGraphCut<float> graphcut;
  graphcut.SetInput(numberofinput,input,numberofatlas+1,datacost);
  graphcut.SetOutput(threshold);
  graphcut.Setdt(timeweight);

  graphcut.SetMode(cutmode);
  graphcut.Run(dataweight,connectmode);

  //finalize
  delete []test;
  delete []tmpatlas;
  delete []imageoffset;
  delete []datacost;
  delete interpolator;
}

void irtkSegmentationFunction::Normalize(irtkRealImage *target, irtkRealImage *source, int n)
{
  double *targetmean, *targetvar, *targetc,
  *sourcemean, *sourcevar, *sourcec,
  *b,*a,*normal;
  irtkRealPixel min, max;
  // Default parameters
  double treshold = 0.0001;
  int iterations = 50;
  int i,k;
  b = new double[n];
  a = new double[n];
  normal = new double[n];
  //Default settings for Gaussian Mixture parameters
  targetmean = new double[n];
  targetvar  = new double[n];
  targetc    = new double[n];

  target->GetMinMax(&min,&max);

  for (i=0;i<n;i++)  targetmean[i] = min + i*(max-min)/(double) n;
  for (i=0;i<n;i++)  targetvar[i] = ((max-min)/(double) n)*((max-min)/(double) n);
  for (i=0;i<n;i++)  targetc[i] = 1.0/(double) n;

  irtkEMClassification classification;
  classification.SetInput(*target);
  classification.SetPadding(-1);
  classification.CreateMask();
  classification.InitialiseGMMParameters(n,targetmean,targetvar,targetc);

  double rel_diff;
  i=1;
  cout << "Segmenting target tissues"<<endl;
  do {
    ///cout << "Iteration = " << i << " / " << iterations << endl;
    rel_diff = classification.IterateGMM(i);
    i++;
  } while ((rel_diff>treshold)&&(i<iterations));
  classification.GetMean(targetmean);
  classification.GetVariance(targetvar);

  sourcemean = new double[n];
  sourcevar  = new double[n];
  sourcec    = new double[n];

  source->GetMinMax(&min,&max);

  for (i=0;i<n;i++)  sourcemean[i] = min + i*(max-min)/(double) n;
  for (i=0;i<n;i++)  sourcevar[i] = ((max-min)/(double) n)*((max-min)/(double) n);
  for (i=0;i<n;i++)  sourcec[i] = 1.0/(double) n;

  classification.SetInput(*source);
  classification.SetPadding(-1);
  classification.CreateMask();
  classification.InitialiseGMMParameters(n,sourcemean,sourcevar,sourcec);

  i=1;
  cout << "Segmenting source tissues"<<endl;
  do {
    ///cout << "Iteration = " << i << " / " << iterations << endl;
    rel_diff = classification.IterateGMM(i);
    i++;
  } while ((rel_diff>treshold)&&(i<iterations));
  classification.GetMean(sourcemean);
  classification.GetVariance(sourcevar);

  for(i=0;i<n;i++) {
    b[i] = targetvar[i]/sourcevar[i];
    a[i] = targetmean[i] - sourcemean[i];
  }
  irtkGaussian *gaussian = new irtkGaussian[n];
  for(i=0;i<n;i++) {
    gaussian[i].Initialise(sourcemean[i],sourcevar[i]);
    normal[i] = gaussian[i].Evaluate(sourcemean[i]);
  }

  irtkRealPixel *ptr = source->GetPointerToVoxels();
  for(k=0;k<source->GetNumberOfVoxels();k++) {
    double value = *ptr;
    double newvalue = 0;
    double weight = 0;
    for(i=0;i<n;i++) {
      newvalue += (gaussian[i].Evaluate(value)/normal[i])*((value - sourcemean[i])*b[i]+targetmean[i]);
      weight += gaussian[i].Evaluate(value)/normal[i];
    }
    newvalue = newvalue/weight;
    *ptr = newvalue;
    ptr++;
  }
  target->GetMinMax(&min,&max);
  source->PutMinMax(min,max);
  // now do the transformation

  delete []targetmean;
  delete []targetvar;
  delete []targetc;
  delete []sourcemean;
  delete []sourcevar;
  delete []sourcec;
  delete []a;
  delete []b;
  delete []normal;
  delete []gaussian;
}

void irtkSegmentationFunction::DetectBackGround(irtkRealImage *target, irtkRealImage *source, int n)
{
  double *targetmean, *targetvar, *targetc,
  *sourcemean, *sourcevar, *sourcec;
  irtkRealPixel min, max;
  // Default parameters
  double treshold = 0.0001;
  int iterations = 50;
  int i,k;
  //Default settings for Gaussian Mixture parameters
  targetmean = new double[n];
  targetvar  = new double[n];
  targetc    = new double[n];

  target->GetMinMax(&min,&max);

  for (i=0;i<n;i++)  targetmean[i] = min + i*(max-min)/(double) n;
  for (i=0;i<n;i++)  targetvar[i] = ((max-min)/(double) n)*((max-min)/(double) n);
  for (i=0;i<n;i++)  targetc[i] = 1.0/(double) n;

  irtkEMClassification classification;
  classification.SetInput(*target);
  classification.SetPadding(-1);
  classification.CreateMask();
  classification.InitialiseGMMParameters(n,targetmean,targetvar,targetc);

  double rel_diff;
  i=1;
  cout << "Segmenting target tissues"<<endl;
  do {
    ///cout << "Iteration = " << i << " / " << iterations << endl;
    rel_diff = classification.IterateGMM(i);
    i++;
  } while ((rel_diff>treshold)&&(i<iterations));
  classification.ConstructSegmentationNoBG(*target);

  sourcemean = new double[n];
  sourcevar  = new double[n];
  sourcec    = new double[n];

  source->GetMinMax(&min,&max);

  for (i=0;i<n;i++)  sourcemean[i] = min + i*(max-min)/(double) n;
  for (i=0;i<n;i++)  sourcevar[i] = ((max-min)/(double) n)*((max-min)/(double) n);
  for (i=0;i<n;i++)  sourcec[i] = 1.0/(double) n;

  irtkEMClassification sclassification;
  sclassification.SetInput(*source);
  sclassification.SetPadding(-1);
  sclassification.CreateMask();
  sclassification.InitialiseGMMParameters(n,sourcemean,sourcevar,sourcec);

  i=1;
  cout << "Segmenting source tissues"<<endl;
  do {
    ///cout << "Iteration = " << i << " / " << iterations << endl;
    rel_diff = sclassification.IterateGMM(i);
    i++;
  } while ((rel_diff>treshold)&&(i<iterations));
  sclassification.ConstructSegmentationNoBG(*source);

  irtkRealPixel *ptr = target->GetPointerToVoxels();
  for(k=0;k<target->GetNumberOfVoxels();k++) {
    if(*ptr == 1 || *ptr == n) {
      *ptr = 1;
    } else {
      *ptr = 0;
    }
    ptr++;
  }

  ptr = source->GetPointerToVoxels();
  for(k=0;k<source->GetNumberOfVoxels();k++) {
    if(*ptr == 1 || *ptr == n) {
      *ptr = 1;
    } else {
      *ptr = 0;
    }
    ptr++;
  }
  // now do the transformation

  delete []targetmean;
  delete []targetvar;
  delete []targetc;
  delete []sourcemean;
  delete []sourcevar;
  delete []sourcec;
}

void irtkSegmentationFunction::EvaluateGraphCut (irtkGreyImage *threshold, irtkGreyImage **input, int number,double timeweight,int cutmode,int connectmode, double dataweight,int padding, int n)
{
  double treshold = 0.0001;
  int iterations = 50;

  irtkEMClassification *classification;
  classification= new irtkEMClassification();
  double *mean, *var, *c;
  double rel_diff;

  //Default settings for Gaussian Mixture parameters
  mean = new double[n];
  var  = new double[n];
  c    = new double[n];

  irtkGreyPixel min, max;
  input[0]->GetMinMax(&min,&max);

  for (int i=0;i<n;i++)  mean[i] = min + i*(max-min)/(double) n;
  for (int i=0;i<n;i++)  var[i] = ((max-min)/(double) n)*((max-min)/(double) n);
  for (int i=0;i<n;i++)  c[i] = 1.0/(double) n;

  classification->SetInput(*input[0]);
  classification->SetPadding(padding);
  classification->CreateMask();
  classification->InitialiseGMMParameters(n,mean,var,c);

  int i=1;
  do {
    cout << "Iteration = " << i << " / " << iterations << endl;
    rel_diff = classification->IterateGMM(i);
    i++;
  } while ((rel_diff>treshold)&&(i<iterations));

  //weights for graph cut
  irtkRealImage **atest;
  atest=new irtkRealImage*[n];
  for (int i=0;i<n;i++){
      atest[i] = new irtkRealImage(input[0]->GetImageAttributes());
      classification->GetProbMap(i,*atest[i]);
  }

  irtkImageGraphCut<irtkGreyPixel> graphcut;
  graphcut.SetInput(number,input,n,atest);
  graphcut.SetOutput(threshold);
  graphcut.SetMode(cutmode);
  graphcut.Setdt(timeweight);
  graphcut.Run(dataweight,connectmode);

  for (int i=0;i<n;i++){
    delete atest[i];
  }
  delete []atest;
  delete []mean;
  delete []var;
  delete []c;
  delete classification;
}

void irtkSegmentationFunction::EvaluateThreshold (irtkRealImage *threshold, irtkGreyImage *target, int padding, int n)
{
  double *mean, *var, *c;
  irtkGreyPixel min, max;
  // Default parameters
  double treshold = 0.0001;
  int iterations = 50;
  int i,j,k;
  //Default settings for Gaussian Mixture parameters
  mean = new double[n];
  var  = new double[n];
  c    = new double[n];

  target->GetMinMax(&min,&max);

  for (i=0;i<n;i++)  mean[i] = min + i*(max-min)/(double) n;
  for (i=0;i<n;i++)  var[i] = ((max-min)/(double) n)*((max-min)/(double) n);
  for (i=0;i<n;i++)  c[i] = 1.0/(double) n;

  irtkEMClassification classification;
  classification.SetInput(*target);
  classification.SetPadding(padding);
  classification.CreateMask();
  classification.InitialiseGMMParameters(n,mean,var,c);

  double rel_diff, tmp_diff;
  i=1; rel_diff = 1;
  cout << "Segmenting tissues: type Background Tissue Blood"<<endl;
  do {
    tmp_diff = rel_diff;
    ///cout << "Iteration = " << i << " / " << iterations << endl;
    rel_diff = classification.IterateGMM(i);
    i++;
  } while ((rel_diff>treshold)&&(i<iterations)&&(tmp_diff > rel_diff));
  cout << "Segmentation done, threshold matrix outputted"<<endl;
  classification.ConstructSegmentationNoBG(*threshold);
  //threshold->Write("d:\\threshold.gipl");
  //remap threshold
  if( n != 3) {
    for(k = 0; k < threshold->GetZ(); k++) {
      for(j = 0; j < threshold->GetY(); j++) {
        for(i = 0; i < threshold->GetX(); i++) {
          if(threshold->GetAsDouble(i,j,k) == 2) {
            threshold->PutAsDouble(i,j,k,3);
          }
          if(threshold->GetAsDouble(i,j,k) == 4) {
            threshold->PutAsDouble(i,j,k,2);
          }
          if(threshold->GetAsDouble(i,j,k) == 5) {
            threshold->PutAsDouble(i,j,k,2);
          }
        }
      }
    }
  } else {
    for(k = 0; k < threshold->GetZ(); k++) {
      for(j = 0; j < threshold->GetY(); j++) {
        for(i = 0; i < threshold->GetX(); i++) {
          if(threshold->GetAsDouble(i,j,k) == 2) {
            threshold->PutAsDouble(i,j,k,3);
          } else if(threshold->GetAsDouble(i,j,k) == 3) {
            threshold->PutAsDouble(i,j,k,2);
          }
        }
      }
    }
  }

  delete mean;
  delete var;
  delete c;
}