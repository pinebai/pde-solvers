function [tcomp, q, qAll] = odeLoad2Buffer( path, filename )

eval(['load ', path,'/', filename, '.dat']);
eval(['qAll = ',filename,';']);
eval(['clear ',filename]);

tcomp = qAll(:,1);
q.c   = qAll(:,2);
q.h   = qAll(:,3);
q.b1  = qAll(:,4);
q.b2  = qAll(:,5);
q.fc  = qAll(:,6);
q.fh  = qAll(:,7);
q.fb1  = qAll(:,8);
q.fb2  = qAll(:,9);

