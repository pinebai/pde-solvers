#include "GenericGraphicsInterface.h"
#include "HDF_DataBase.h"
#include "ShowFileReader.h"
#include "DisplayParameters.h"
#include "display.h"
#include "DerivedFunctions.h"
#include "fileOutput.h"
#include "PlotIt.h"

static int totalNumberOfArrays=0;

void
checkArrays(const aString & label)
// Output a warning messages if the number of arrays has increased
{
  if( false )  // **** set this to true to turn on array checking ***
  {
    if(GET_NUMBER_OF_ARRAYS > totalNumberOfArrays ) 
    {
      totalNumberOfArrays=GET_NUMBER_OF_ARRAYS;
      printf("\n**** %s :Number of A++ arrays = %i \n\n",(const char*)label,GET_NUMBER_OF_ARRAYS);
    }
  }
}

static int 
buildMainMenu( aString *menu0,
               aString *&menu,
               RealGridCollectionFunction & u,
               aString *sequenceName,
               const int & numberOfSolutions,
               const int & numberOfComponents,
               const int & numberOfSequences,
               int & chooseAComponentMenuItem,
               int & chooseASolutionMenuItem,
               int & numberOfSolutionMenuItems,
               int & chooseASequenceMenuItem,
               int & numberOfSequenceMenuItems,
               const int & maxMenuSolutions,
	       const int & maximumNumberOfSolutionsInTheMenu,
	       int & solutionIncrement,
	       const int & maxMenuSequences,
	       const int & maximumNumberOfSequencesInTheMenu,
	       int & sequenceIncrement )
{
  // create the real menu by adding in the component names, these will appear as a 
  // cascaded menu
  char buff[120];

  chooseAComponentMenuItem=0;  // menu[chooseAComponentMenuItem]=">choose a component"
  chooseASolutionMenuItem=0;  
  numberOfSolutionMenuItems=0;
  chooseASequenceMenuItem=0;  
  numberOfSequenceMenuItems=0;

  int numberOfMenuItems0=0;
  while( menu0[numberOfMenuItems0]!="" )
  {
    numberOfMenuItems0++;
  }  
  numberOfMenuItems0++;

  // const int maxMenuSolutions=25;  // cascade solution menu if there are more than this many solutions
  // const int maximumNumberOfSolutionsInTheMenu=400;  // stride through the solutions if there are more
  // than this many solutions.
  // int solutionIncrement=1;                          // Here is the stride.

  // const int maxMenuSequences=25;  // cascade sequence menu if there are more than this many sequences
  // const int maximumNumberOfSequencesInTheMenu=400;  // stride through the sequences if there are more
  // than this many sequences.
  // int sequenceIncrement=1;                          // Here is the stride.

  delete [] menu;
  menu = new aString [numberOfMenuItems0+numberOfComponents+numberOfSolutions
		    +numberOfSolutions/maxMenuSolutions 
		    +numberOfSequences+numberOfSequences/maxMenuSequences +2];

  int i=-1;
  for( int i0=0; i0<numberOfMenuItems0 ; i0++ )
  {
    menu[++i]=menu0[i0];    
    if( menu[i]==">choose a component" )
    {
      chooseAComponentMenuItem=i;
      for( int j=0; j<numberOfComponents; j++ )
      {
	menu[++i]=u.getName(u.getComponentBase(0)+j);
	if( menu[i] == "" || menu[i]==" " )
	  menu[i]=sPrintF(buff,"component%i",u.getComponentBase(0)+j);
      }
    }
    else if( menu[i]=="<>choose a solution" )
    {
      // make menu items that display all the solutions. If there many solutions then we cascade the solutions
      // into groups, each group has maxMenuSolutions entries
      chooseASolutionMenuItem=i;
      solutionIncrement=1;
      if( numberOfSolutions>maximumNumberOfSolutionsInTheMenu )
	solutionIncrement=(numberOfSolutions+maximumNumberOfSolutionsInTheMenu-1)/maximumNumberOfSolutionsInTheMenu;
	  
      int k=0;
      for( int j=0; j<numberOfSolutions; j+=solutionIncrement )
      {
	if( numberOfSolutions>maxMenuSolutions && ( k % maxMenuSolutions==0) )
	{
	  if( j==0 )
	    menu[++i]=sPrintF(buff,">solutions %i to %i",j,j+maxMenuSolutions*solutionIncrement-1);
	  else
	    menu[++i]=sPrintF(buff,"<>solutions %i to %i",j,
			      min(j+maxMenuSolutions*solutionIncrement-1,numberOfSolutions-1));
	}
	menu[++i]=sPrintF(buff,"solution%i",j);
	k++;
      }
      if( numberOfSolutions>maxMenuSolutions )
	menu[++i]="< ";
      numberOfSolutionMenuItems=i-chooseASolutionMenuItem;
    }
    else if( menu[i]==">sequence" )
    {
      // make menu items that display all the sequences. If there many sequences then we cascade them into groups.
      chooseASequenceMenuItem=i;
      sequenceIncrement=1;
      if( numberOfSequences>maximumNumberOfSequencesInTheMenu )
	sequenceIncrement=(numberOfSequences+maximumNumberOfSequencesInTheMenu-1)/maximumNumberOfSequencesInTheMenu;
	  
      int k=0;
      for( int j=0; j<numberOfSequences; j+=sequenceIncrement )
      {
	if( numberOfSequences>maxMenuSequences && ( k % maxMenuSequences==0) )
	{
	  if( j==0 )
	    menu[++i]=sPrintF(buff,">sequences %i to %i",j,j+maxMenuSequences*sequenceIncrement-1);
	  else
	    menu[++i]=sPrintF(buff,"<>sequences %i to %i",j,
			      min(j+maxMenuSequences*sequenceIncrement-1,numberOfSequences-1));
	}
	menu[++i]=sequenceName[j];
	k++;
      }
      if( numberOfSequences>maxMenuSequences )
	menu[++i]="< ";
      numberOfSequenceMenuItems=i-chooseASequenceMenuItem;
    }
  }
  return 0;
}




//===============================================================================================
// This plotCellWave program is used to display CellWave simulations stored in a show file
//
//===============================================================================================
//\begin{>plotCellWaveMenuInclude.tex}{}
//\no function header:
//
// Here is a desciption of the menu options available for plotCellWave.
//\begin{description}\index{show files!plotting}
//  \item[contour] : plot contours.
//  \item[stream lines] : plot streamlines.
//  \item[grid] : plot the grid.
//  \item[sequence] : plot a sequence.
//  \item[next] : plot solutions from the next frame.
//  \item[previous] : plot solutions from the previous frame.
//  \item[choose a component] : plot a different component.
//  \item[choose a solution] : choose a different frame.
//  \item[next component] : plot the next component.
//  \item[previous component] : plot the previous component.
//  \item[derived types] : build new components as functions of the old ones, such as the vorticity from the velocity.
//     One derived types have been created they will appear in the plotStuff component menus.
//  \item[movie] : \index{movie} plot a number of frames in a row.
//  \item[movie and save] : plot frames and save each one as a hard copy.
//  \item[plot bounds] : \index{plot bounds}change the manner in which the plot bounds are determined.
//    \begin{description}
//    \item[set plot bounds] : specify bounds for plotting.
//    \item[use default plot bounds] : use default plot bounds.
// \end{description}
//  \item[check mappings with grid] : a debugging option, check validity of mappings after they have been read
//     in from a data-base file.
//  \item[erase] : erase anything in the window.
//  \item[redraw] : redraw the screen.
//  \item[open a new file] : open a new show file for reading.
//  \item[file output] : \index{file output} output results to a file.
//  \item[help] : print a short help list.
//  \item[exit] : exit this menu and continue on (same as 'continue').
// \end{description}
//
//
//\end{plotCellWaveMenuInclude.tex} 


int 
plotCellWave(int argc, char *argv[])
{
  Overture::start(argc,argv);  // initialize Overture

  cout << "Type: `plotCellWave fileName [file.cmd]' to read the show file called fileName \n";
  cout << "       and optionally start reading a command file \n";

  aString nameOfShowFile, commandFileName;
  nameOfShowFile="";
  if( argc > 1 )
  {
    nameOfShowFile=argv[1];
    int l=nameOfShowFile.length()-1;
    if( l>2 && nameOfShowFile(l-3,l)==".cmd" )
    {
      commandFileName=nameOfShowFile;
      nameOfShowFile="";
    }
    else if( argc>2 )
      commandFileName=argv[2];
  }
  
  bool done=FALSE;
    
// create a Graphics Interface
  GenericGraphicsInterface & ps 
    = *Overture::getGraphicsInterface("CellWave visualization",TRUE); 
  GraphicsParameters psp;       // create an object that is used to pass parameters
  aString movieFileName;
    
  // By default start saving the command file called "plotStuff.cmd"
  aString logFile="out_plotCellWave.cmd";
  ps.saveCommandFile(logFile);
  cout << "User commands are being saved in the file `" << (const char *)logFile << "'\n";

  // read from a command file if given
  if( commandFileName!="" )
    ps.readCommandFile(commandFileName);

  checkArrays("plotStuff: before loop");
  aString *sequenceName = NULL;

  // this loop is used to look at more than one file
  while( !done )
  {
    // cout << ">> Enter the name of the show file:" << endl;
    // cin >> nameOfShowFile;
    if( nameOfShowFile=="" )
      ps.inputString(nameOfShowFile,">> Enter the name of the show file:");

    ShowFileReader showFileReader(nameOfShowFile);

    int numberOfFrames=showFileReader.getNumberOfFrames();
    int numberOfSolutions = max(1,numberOfFrames);
    int numberOfComponents;
    int numberOfSequences=showFileReader.getNumberOfSequences();
    if( numberOfSequences>0 )
    {
      delete [] sequenceName;
      sequenceName = new aString[numberOfSequences];
      showFileReader.getSequenceNames(sequenceName,numberOfSequences);
    }
    

    CompositeGrid cg;  

    // set up a function for contour plotting:
    Range all;
    realCompositeGridFunction u;

    char buff[120];
    aString answer,answer2;
    aString menu0[]= { "!plotStuff",
		       "contour",
		       "line plot",
		       "probe plot",
		       "stream lines",
		       "grid",
		       ">sequence",
		       "<next",
		       "previous",
		       ">choose a component", 
		       "<>choose a solution", 
		       "<next component",
		       "previous component", 
		       "derived types",
		       "movie",
		       "movie and save",
		       ">plot bounds",
		       "set plot bounds",
		       "use default plot bounds",
		       "<check mappings with grid",
		       "erase",
		       "redraw",
		       "open a new file",
		       "file output",
		       "help",
		       "exit",
		       "" };
    aString help[]= { 
       "contour                    : plot contours (surfaces)",
       "line plot                  : plot line cut across solution field",
       "probe plot                 : plot timehistory of field at probe locations",
       "stream lines               : draw stream lines",
       "grid                       : plot the grid",
       "sequence                   : plot a sequence that has been saved in the show file",
       "next                       : plot the next solution of all items on the screen",
       "previous                   : plot the previous solution of all items on the screen" ,
       "choose a component         : plot a different component of all items on the screen", 
       "choose a solution          : plot a different solution of all items on the screen", 
       "next component             : plot the next component of all items on the screen",
       "previous component         : plot the previous component of all items on the screen", 
       "derived types              : define new quantities such as vorticity, derivatives etc.",
       "movie                      : plot the next `n' solutions",
       "movie and save             : plot the next `n' solutions and save each as a postscript file",
       "set plot bounds            : specify fixed bounds for plotting. Useful for movies.",
       "use default plot bounds    : let plotStuff determine the plotting bounds",
       "check mappings with grid   : call the checkMapping routine",
       "erase                      : erase everything",
       "redraw                     : force a redraw (useful to call from command files)",
       "open a new file            : open a new show file to read",
       "file output                : output solutions to a file (ascii)",
       "help",
       "exit",
       "" };

    int plotOptions = 0;
    bool plotNewFunction = FALSE;
    bool plotNewComponent= FALSE;
    bool movieMode=FALSE;
    int numberOfMovieFrames=numberOfFrames;
    int solutionNumber=1;
    int component = 0;
    int numberOfHeaderComments;
    
    const aString *headerComment; // This array holds the comments that should go in the plot
    showFileReader.getASolution(solutionNumber,cg,u);
    headerComment=showFileReader.getHeaderComments(numberOfHeaderComments);
    numberOfComponents=u.getComponentDimension(0);
    const int numberOfComponents0=numberOfComponents;
    
    // this next class knows how to form derived quantities such as vorticity, derivatives etc.
    DerivedFunctions derivedFunctions(showFileReader);

    // create the real menu by adding in the component names, these will appear as a 
    // cascaded menu

    int chooseAComponentMenuItem;  // menu[chooseAComponentMenuItem]=">choose a component"
    int chooseASolutionMenuItem;  
    int numberOfSolutionMenuItems=0;
    int chooseASequenceMenuItem;  
    int numberOfSequenceMenuItems=0;
    const int maxMenuSolutions=25;  // cascade solution menu if there are more than this many solutions
    const int maximumNumberOfSolutionsInTheMenu=400;  // stride through the solutions if there are more
                                                      // than this many solutions.
    int solutionIncrement=1;                          // Here is the stride.

    const int maxMenuSequences=25;  // cascade sequence menu if there are more than this many sequences
    const int maximumNumberOfSequencesInTheMenu=400;  // stride through the sequences if there are more
                                                      // than this many sequences.
    int sequenceIncrement=1;                          // Here is the stride.

    aString *menu=NULL;
    buildMainMenu( menu0,
		   menu,
                   u,
                   sequenceName,
                   numberOfSolutions,
                   numberOfComponents,
                   numberOfSequences,
		   chooseAComponentMenuItem,
		   chooseASolutionMenuItem,
		   numberOfSolutionMenuItems,
		   chooseASequenceMenuItem,
		   numberOfSequenceMenuItems,
		   maxMenuSolutions,
		   maximumNumberOfSolutionsInTheMenu,
                   solutionIncrement,
		   maxMenuSequences,
		   maximumNumberOfSequencesInTheMenu,
		   sequenceIncrement );


    psp.set(GI_TOP_LABEL,headerComment[0]);  // set title
    psp.set(GI_TOP_LABEL_SUB_1,headerComment[1]);  
    psp.set(GI_TOP_LABEL_SUB_2,headerComment[2]);  
    if( cg.numberOfDimensions()==1 )
      psp.set(GI_COLOUR_LINE_CONTOURS,TRUE);
    
    // set default prompt
    ps.appendToTheDefaultPrompt("plotCellWave>");

    int menuItem=-1;
    for( int it=0; ; it++)
    {
      checkArrays("plotCellWave: in for(;;)");

      if( it==0 && numberOfFrames<=0 )
        answer="grid";
      else
        menuItem=ps.getMenuItem(menu,answer);
      if( answer=="grid" )
      {
	PlotIt::plot(ps, cg, psp);   // plot the composite grid
	if( psp.getObjectWasPlotted() ) 
	  plotOptions |= 1;

	if( false )
	{
	  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
	    cg[grid].displayComputedGeometry();
	}
      }
      else if( answer=="contour" )
      {
	PlotIt::contour(ps, u, psp);  // contour/surface plots
	if( psp.getObjectWasPlotted() & 1 ) 
	  plotOptions |= 2;
	if( psp.getObjectWasPlotted() & 2 )  // grid was also plotted
	  plotOptions |= 1;
      }
      else if( answer=="line plot" )
      {
	printf(" -- ERROR: line plot not available yet. --\n");

	int oldPBGG = ps.getPlotTheBackgroundGrid();
	int oldKAR  = ps.getKeepAspectRatio();
	const GridCollection & gc = *(u.gridCollection);
	RealArray xBound(2,3);
	PlotIt::getPlotBounds(gc,psp,xBound);
	
	// plot solution on lines that cut the 2D grid
	//contourCuts(gi, uGCF,psp );
	PlotIt::contourCuts(ps, u,  psp );
	
	// Restore plotbackgroundgrid and keepAspectRatio after this call
	ps.setPlotTheBackgroundGrid(oldPBGG);
	
	psp.keepAspectRatio=oldKAR;
	ps.setKeepAspectRatio(psp.keepAspectRatio); 

	// the boundingbox is messed up (set for 1D) after this call
	ps.setGlobalBound(xBound);
	
	// erase the labels and replot them
	ps.eraseLabels(psp);

	// replot the 3D object
	//plotObject = TRUE;
	//plotContours = TRUE;

	//PlotIt::contour(ps, u, psp);  // contour/surface plots
	//if( psp.getObjectWasPlotted() & 1 ) 
	//  plotOptions |= 2;
	//if( psp.getObjectWasPlotted() & 2 )  // grid was also plotted
	//  plotOptions |= 1;
      }
      else if( answer=="probe plot" )
      {
	printf(" -- ERROR: probe plot not available yet. --\n");
	//	PlotIt::contour(ps, u, psp);  // contour/surface plots
	//if( psp.getObjectWasPlotted() & 1 ) 
	//  plotOptions |= 2;
	//if( psp.getObjectWasPlotted() & 2 )  // grid was also plotted
	//  plotOptions |= 1;
      }
      else if( answer=="stream lines" )
      {
	PlotIt::streamLines(ps, u, psp);  // streamlines
	if( psp.getObjectWasPlotted() ) 
	  plotOptions |= 4;
      }
      else if( answer=="derived types" )
      {
        if( numberOfComponents>0 )
	{
	  aString *componentNames = new aString [numberOfComponents];
	  for( int n=0; n<numberOfComponents; n++ )
	    componentNames[n]=u.getName(n);
	  
	  derivedFunctions.update(ps,numberOfComponents,componentNames);
          delete [] componentNames;

          derivedFunctions.getASolution(solutionNumber,cg,u);

          numberOfComponents=numberOfComponents0+derivedFunctions.numberOfDerivedTypes();
	  buildMainMenu( menu0,
			 menu,
			 u,
			 sequenceName,
			 numberOfSolutions,
			 numberOfComponents, 
			 numberOfSequences,
			 chooseAComponentMenuItem,
			 chooseASolutionMenuItem,
			 numberOfSolutionMenuItems,
			 chooseASequenceMenuItem,
			 numberOfSequenceMenuItems,
			 maxMenuSolutions,
			 maximumNumberOfSolutionsInTheMenu,
			 solutionIncrement,
			 maxMenuSequences,
			 maximumNumberOfSequencesInTheMenu,
			 sequenceIncrement );


	}
	else
	{
	  printf("ERROR: no components are available\n");
	}
      }
      else if( menuItem > chooseASequenceMenuItem && menuItem <= chooseASequenceMenuItem+numberOfSequenceMenuItems )
      {
        // plot a sequence
	int sequenceNumber=menuItem-chooseASequenceMenuItem;
        if( numberOfSequences>maxMenuSequences )
	{
          // adjust the sequence number when there are many sequences since we add in extra menu items
          // into the list.
          int extra = 1+ sequenceNumber/maxMenuSequences;  
          extra=1+ (sequenceNumber-extra)/maxMenuSequences;
          sequenceNumber-=extra;
          // printf("menuItem-chooseASequenceMenuItem=%i, extra=%i sequenceNumber=%i\n",
          //     menuItem-chooseASequenceMenuItem,extra,sequenceNumber);
	}
        sequenceNumber=(sequenceNumber-1)*sequenceIncrement;

        assert( sequenceNumber>=0 && sequenceNumber<numberOfSequences );
	aString name;
	realArray time,value;
	const int maxComponentName1=25, maxComponentName2=1;
	aString componentName1[maxComponentName1], componentName2[maxComponentName2];
	  
	showFileReader.getSequence(sequenceNumber,name,time,value,
				   componentName1,maxComponentName1,
				   componentName2,maxComponentName2);
	// printf("sequence %i: name=%s\n",sequenceNumber,(const char*)name);
        // display(value,"value");
	
	  
	ps.erase();
	psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,FALSE);
	psp.set(GI_TOP_LABEL_SUB_1,"");
	psp.set(GI_TOP_LABEL_SUB_2,"");
        Range all;
	PlotIt::plot(ps, time, value(all,all,value.getBase(2)), name, "t", componentName1, psp);
	psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,TRUE);
	ps.erase();

      }
      else if( answer=="next" )
      {
	solutionNumber = (solutionNumber % numberOfSolutions) +1;
	plotNewFunction=TRUE;
      }
      else if( answer=="previous" )
      {
	solutionNumber = ((solutionNumber-2+numberOfSolutions) % numberOfSolutions) +1;
	plotNewFunction=TRUE;
      }
      else if( answer=="next component" )
      {
	component= (component+1) % numberOfComponents;
	plotNewComponent=TRUE;
      }
      else if( answer=="previous component" )
      {
	component= (component-1+numberOfComponents) % numberOfComponents;
	plotNewComponent=TRUE;
      }
      else if( menuItem > chooseAComponentMenuItem && menuItem <= chooseAComponentMenuItem+numberOfComponents )
      {
	component=menuItem-chooseAComponentMenuItem-1 + u.getComponentBase(0);
	plotNewComponent=TRUE;
	// cout << "chose component number=" << component << endl;
      }
      else if( answer=="choose a component" )
      { // *** not used *** Make a menu with the component names. If there are no names then use the component numbers
	aString *menu2 = new aString[numberOfComponents+1];
	for( int i=0; i<numberOfComponents; i++ )
	{
	  menu2[i]=u.getName(u.getComponentBase(0)+i);
	  if( menu2[i] == "" || menu2[i]==" " )
	    menu2[i]=sPrintF(buff,"component%i",u.getComponentBase(0)+i);
	}
	menu2[numberOfComponents]="";   // null string terminates the menu
	component = ps.getMenuItem(menu2,answer2);
	component+=u.getComponentBase(0);
	delete [] menu2;
	plotNewComponent=TRUE;
      }
      else if( menuItem > chooseASolutionMenuItem && menuItem <= chooseASolutionMenuItem+numberOfSolutionMenuItems )
      {
	solutionNumber=menuItem-chooseASolutionMenuItem;
        if( numberOfSolutions>maxMenuSolutions )
	{
          // adjust the solution number when there are many solutions since we add in extra menu items
          // into the list.
          int extra = 1+ solutionNumber/maxMenuSolutions;  
          extra=1+ (solutionNumber-extra)/maxMenuSolutions;
          solutionNumber-=extra;
          // printf("menuItem-chooseASolutionMenuItem=%i, extra=%i solutionNumber=%i\n",
          //     menuItem-chooseASolutionMenuItem,extra,solutionNumber);
	}
        solutionNumber=(solutionNumber-1)*solutionIncrement+1;
	
	plotNewFunction=TRUE;
      }
      else if( answer=="choose a solution" )
      { // ***** not used ****  Make a menu with the solution Names
	aString *menu2 = new aString[numberOfFrames+1];
	for( int i=0; i<numberOfFrames; i++ )
	  menu2[i]=sPrintF(buff,"solution%i",i);
	menu2[numberOfFrames]="";   // null string terminates the menu
	solutionNumber = ps.getMenuItem(menu2,answer2)+1;
	delete [] menu2;
	plotNewFunction=TRUE;
      }
      else if( answer=="movie" || answer=="movie and save" )
      {
	movieMode=TRUE;
	numberOfMovieFrames=numberOfFrames;
	ps.inputString(answer2,sPrintF(buff,"Enter the number of frames (total=%i)",numberOfFrames));
	if( answer2 !="" && answer2!=" ")
	{
	  sScanF(answer2,"%i",&numberOfMovieFrames);
	  printf("number of frames = %i \n",numberOfMovieFrames);
	}
        if( answer=="movie and save" )
	{
	  ps.inputString(answer2,"Enter basic name for the ppm files (default=plot)");
	  if( answer2 !="" && answer2!=" ")
	    movieFileName=answer2;
          else
	    movieFileName="plot";
          ps.outputString(sPrintF(buff,"pictures will be named %s0.ppm, %s1.ppm, ...",
            (const char*)movieFileName,(const char*)movieFileName));
	}
	
      }
      else if( answer=="set plot bounds" )
      {
        RealArray xBound(2,3);
	xBound=0.;
        xBound(1,Range(0,2))=1.;
        if( cg.numberOfDimensions()==2 )
  	  ps.inputString(answer2,sPrintF(buff,"Enter bounds xa,xb, ya,yb "));
        else
  	  ps.inputString(answer2,sPrintF(buff,"Enter bounds xa,xb, ya,yb, za,zb "));
        if( answer2!="" )
 	  sScanF(answer2,"%e %e %e %e %e %e",&xBound(0,0),&xBound(1,0),&xBound(0,1),&xBound(1,1),
              &xBound(0,2),&xBound(1,2));
	
        ps.resetGlobalBound(ps.getCurrentWindow());
	ps.setGlobalBound(xBound);
	
	psp.set(GI_PLOT_BOUNDS,xBound); // set plot bounds
	psp.set(GI_USE_PLOT_BOUNDS,TRUE);  // use the region defined by the plot bounds
      }
      else if( answer=="use default plot bounds" )
      {
        psp.set(GI_USE_PLOT_BOUNDS,FALSE);  // use the region defined by the plot bounds
      }
      else if( answer=="check mappings with grid" )
      {
	for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
	{
	  if( cg[grid].mapping().mapPointer==NULL )
	  {
	    cout << "ERROR: This grid has no mappings! \n";
	    break;
	  }
	  cg[grid].mapping().checkMapping();
	}
      }
      else if( answer=="erase" )
      {
	ps.erase();
	plotOptions=0;
      }
      else if( answer=="redraw" )
      { // force a redraw -- add to command files to force the drawing of the screen
	ps.redraw(TRUE);
      }
      else if( answer=="open a new file" )
      {
        nameOfShowFile=""; // do this so we prompt for a new name
	break;
      }
      else if( answer=="file output" )
      {
        fileOutput(ps, u);
      }
      else if( answer=="exit" )
      {
        done=TRUE;
	break;
      }
      else if( answer=="help" )
      {
	for( int i=0; help[i]!=""; i++ )
	   ps.outputString(help[i]);
      }
      else
      {
        cout << "unknown response, answer=[" << answer << "]\n";
	ps.stopReadingCommandFile();
      }
      if( movieMode )
      { // ************** Movie Mode *******************
	psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,TRUE);
	for( int frame=1; frame<=numberOfMovieFrames; frame++ )
	{
          if( answer=="movie and save" )
	  { // save a ppm file
            psp.set(GI_HARD_COPY_TYPE,GraphicsParameters::ppm);
            ps.outputString(sPrintF(buff,"Saving file %s%i.ppm",(const char*)movieFileName,frame-1));
	    ps.hardCopy(    sPrintF(buff,            "%s%i.ppm",(const char*)movieFileName,frame-1),psp);
            psp.set(GI_HARD_COPY_TYPE,GraphicsParameters::postScript);
	  }

	  solutionNumber = (solutionNumber+numberOfSolutions) % numberOfSolutions +1;

          // showFileReader.getASolution(solutionNumber,cg,u);
          derivedFunctions.getASolution(solutionNumber,cg,u);

          headerComment=showFileReader.getHeaderComments(numberOfHeaderComments);
          numberOfComponents=u.getComponentDimension(0);

	  psp.set(GI_TOP_LABEL,headerComment[0]);  // set title
	  psp.set(GI_TOP_LABEL_SUB_1,headerComment[1]);  
	  psp.set(GI_TOP_LABEL_SUB_2,headerComment[2]);  
	  ps.erase();
	  if( plotOptions & 1 )
	    PlotIt::plot(ps, cg, psp );
	  if( plotOptions & 2 )
	    PlotIt::contour(ps, u, psp );
	  if( plotOptions & 4 )
	    PlotIt::streamLines(ps, u, psp ); 

	  ps.redraw(TRUE);   // *****
	}
	if( answer=="movie and save" )
	{ // save a ppm file
          psp.set(GI_HARD_COPY_TYPE,GraphicsParameters::ppm);
	  ps.outputString(sPrintF(buff,"Saving file %s%i.ppm",(const char*)movieFileName,numberOfMovieFrames));
	  ps.hardCopy(    sPrintF(buff,            "%s%i.ppm",(const char*)movieFileName,numberOfMovieFrames),psp);
	  psp.set(GI_HARD_COPY_TYPE,GraphicsParameters::postScript);
	}

	psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,FALSE);
	movieMode=FALSE;
      }      
      else if( plotNewFunction || plotNewComponent )
      {
	if( plotNewFunction )
	{
          // showFileReader.getASolution(solutionNumber,cg,u);
          derivedFunctions.getASolution(solutionNumber,cg,u);

          headerComment=showFileReader.getHeaderComments(numberOfHeaderComments);
          numberOfComponents=u.getComponentDimension(0);

	  psp.set(GI_TOP_LABEL,headerComment[0]);  // set title
	  psp.set(GI_TOP_LABEL_SUB_1,headerComment[1]);  
	  psp.set(GI_TOP_LABEL_SUB_2,headerComment[2]);  
	  plotNewFunction=FALSE;
	}
	if( plotNewComponent )
	{
	  psp.set(GI_COMPONENT_FOR_CONTOURS,component);
	  plotNewComponent=FALSE;
	}
	psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,TRUE);
	ps.erase();
	if( plotOptions & 1 )
	  PlotIt::plot(ps, cg, psp );
	if( plotOptions & 2 )
	  PlotIt::contour(ps, u, psp );
	if( plotOptions & 4 )
	  PlotIt::streamLines(ps, u, psp );

	psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,FALSE);
      }
    } // end for(;;)

    delete [] menu;
    ps.unAppendTheDefaultPrompt(); // reset defaultPrompt

    if( true )
    {
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
	cg[grid].displayComputedGeometry();
    }
  } // end while not done
  
  
  delete [] sequenceName;


  Overture::finish();          
  return 0;
}
