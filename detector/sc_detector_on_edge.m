function [Hypo_Prob] = ...
    sc_detector_on_edge(Obj_Msk,I_edge,codebook,para,verbose)
%SC_DETECTOR_ON_EDGE(I_edge,codebook,para,verbose) pedestrian detection
% using shape context voting
%
% INPUT:
%   Obj_Msk:    Foreground mask
%   I_EDGE:     edge map for different scales
%   CODEBOOK:   object models information
%   PARA:       parameters, please refer to 'set_parameter.m'
%   VERBOSE:    debug level control
% OUTPUT:
%   Hypo_Prob:  Hyposthesis for the person's contours

  global Img
  global PN

%%%%%%%%%%%%%% parameter setup %%%%%%%%%%%%%%%
para_sc = para{1};
para_fea= para{2};
para_vote=para{3};
voter_filter    = para_vote.voter_filter;

win_thresh  = [];
recog_result= [];


%%%%%%%%%%%%%% iterate on scale %%%%%%%%%%%%%%

ratio       = para_fea.ratio;
if(verbose);
    begin_time  = cputime;
    last_time   = begin_time;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % %  basic voting % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
tic

   
    edge_map     = I_edge.edge;
    theta_map    = I_edge.theta;
    recog_res	= sc_vote_1_scale(Obj_Msk,edge_map, theta_map, codebook, para, verbose);
    recog_result= [recog_result,recog_res];
    win_thresh  = [win_thresh; recog_res.win_thresh];


if(ratio>1)
    ratio=1/ratio;
end

max_score_thresh  =   max(win_thresh);
ratio_list  = 1;

% % % % % % % % % % generate hypo mask and bbox,prune by height % % % % % %
recog_result    = generate_mask_and_bbox(recog_result, codebook, ...
    max_score_thresh,para_vote, verbose);



Prob_Carried_Object=recog_result.edge-(recog_result.edge.*recog_result.hypo_mask>0.7);
% % % % % % % % % % % % Edge Link (Find Connected Componet) % % % % %  % %

edgeim=Prob_Carried_Object>0.05;% % % %Previous was 0.15

% Fill Gaps 3 pixels
gapsize=2;

edgeim = filledgegaps(edgeim, gapsize);

minlength=10;

[edgelist] = edgelink(edgeim, minlength);


% Display the edgelists with random colours for each distinct edge
% in figure 2
% % % % % % % % % Find the curve % % % % % % % % %
Nedge = length(edgelist);

TNedge=1;
for L=1:Nedge
    
    Temp_edgelist{TNedge}=edgelist{L};
     y = edgelist{L}(:,1);   % Note that (col, row) corresponds to (x,y)
     x = edgelist{L}(:,2);
     Npts = length(x);
       D = sqrt((x(1)-x(Npts))^2 + (y(1)-y(Npts))^2);
            
            y1my2 = y(1)-y(Npts);  % pluse one to ab                     % Pre-compute parameters
            x2mx1 = x(Npts)-x(1);
            if y1my2==0 | x2mx1==0
                Npts=Npts-1;
                y1my2 = y(1)-y(Npts);  % pluse one to ab                     % Pre-compute parameters
                x2mx1 = x(Npts)-x(1);
            end
            m=y1my2/x2mx1;
            C = y(Npts)*x(1) - y(1)*x(Npts);
            d = (x*y1my2 + y*x2mx1 + C)/D;
            
               if Npts>70 %Npts=70 for adding turn point
                    %%----Smooth the curve -------%%
                     windowSize = 10;
                      b = (1/windowSize)*ones(1,windowSize);
                      a=1;
                      y = filter(b,a,d);
                      %---Find the Exterma----%
                      [ymax,imax,ymin,imin] = extrema(y);
                      En_p=find(imin==1|imin==Npts);
                      imin(En_p)=[];
                      ymin(En_p)=[];
                      Exa=[imax ymax;imin ymin];
                      Exa=sort(Exa);
                      TP=[];
                      
                      for ex=1:size(Exa,1)-1
                          TP=[TP round((Exa(ex,1)+Exa(ex+1,1))/2)];
                      end
                      
                      TP=[1 TP Npts];
                else
                    TP=[1 Npts];
                   % Select the TP in the middle of contour
                  
               end
               if ~isempty(TP)
                NTP=size(TP,2)-1;
                for t1=1:NTP
                     Temp_edgelist{TNedge}=edgelist{L}(TP(t1):TP(t1+1),:);
                     TNedge=TNedge+1;
                    
                end
               else
                   TNedge=TNedge+1;
               end
      
end
edgelist=Temp_edgelist;
Nedge = length(edgelist);
Ne=1;
while Ne<=Nedge
      x = edgelist{Ne}(:,2);
     if size(x,1)<5
         edgelist(Ne)=[];
         Nedge=Nedge-1;
         Ne=Ne-1;
     end
     Ne=Ne+1;
end
Nedge = length(edgelist);
CurveEdge=0;
% remove small edge

for e = 1:Nedge
    
    y = edgelist{e}(:,1);   % Note that (col, row) corresponds to (x,y)
    x = edgelist{e}(:,2);
    [Curve_m] = maxlinedev(x,y);
    
    Curve_List{e}=Curve_m;
    
    Curve_T=0.9; % % %was 0.9 for Pets 
    
    if Curve_m>Curve_T && length(edgelist{e}(:,1))>=15
        if isinf(Curve_m)
           CurveEdge=CurveEdge+1;
            if length(x)>20
            N_Seed=20;
            else
                N_Seed=length(x);
            end
            idx=randperm(length(x));
            Seed1=x(idx(1:N_Seed));
            Seed2=y(idx(1:N_Seed));
            Seed_point_6{CurveEdge}=[Seed1 Seed2];
            New_edgelist{CurveEdge}=edgelist{e};
        else
            CurveEdge=CurveEdge+1;
            New_edgelist{CurveEdge}=edgelist{e};
            Npts = length(x);
            D = sqrt((x(1)-x(Npts))^2 + (y(1)-y(Npts))^2);
            
            y1my2 = y(1)-y(Npts);  % pluse one to ab                     % Pre-compute parameters
            x2mx1 = x(Npts)-x(1);
            if y1my2==0 | x2mx1==0
                Npts=Npts-1;
                y1my2 = y(1)-y(Npts);  % pluse one to ab                     % Pre-compute parameters
                x2mx1 = x(Npts)-x(1);
            end
            m=y1my2/x2mx1;
            C = y(Npts)*x(1) - y(1)*x(Npts);
            d = (x*y1my2 + y*x2mx1 + C)/D;
            idx =(d==0);
            d=(1-idx).*d+idx;
           % if  (m<=1)  &   (m>-4)
                yy1=[];
                yy2=[];
                Seed_point1=[];
                Seed_point2=[];
                % -------Find turn point ----------
                for i=1:Npts
               
                   %      if  (m<=1)  &   (m>-4)                  %%length(x)-length(unique(x))<15 %%%%%%%%%%% Previously was 15..we should solve when we need the other part
                    if abs(d(i,1))>1
                        Bet_d=sign(d(i,1))*[1:round(abs(d(i,1)))];
                        Ld=length(Bet_d);
                        X1=repmat(x(i),[1,Ld]);
                        Y1=round((D*Bet_d-X1*y1my2-C)/x2mx1);
                        yy1=[X1',Y1';yy1];
                        Pnt_Bet=[X1',Y1'];
                        Weight_Pnt=repmat(Pnt_Bet,[round(abs(d(i,1))),1]);
                        Seed_point1=[Seed_point1;Weight_Pnt];
                        
                    else
                        y1=round((D*d(i,1)-x(i)*y1my2-C)/x2mx1);
                        yy1=[x(i),y1;yy1];
                        Seed_point1=[Seed_point1;[x(i),y1]];
                        
                    end
                    %  else
                    
                    if abs(d(i,1))>1
                        Bet_d=sign(d(i,1))*[1:round(abs(d(i,1)))];
                        Ld=length(Bet_d);
                        Y1=repmat(y(i),[1,Ld]);
                        X1=round((D*Bet_d-Y1*x2mx1-C)/y1my2);
                        yy2=[X1',Y1';yy2];
                        Pnt_Bet=[X1',Y1'];
                        Weight_Pnt=repmat(Pnt_Bet,[round(abs(d(i,1))),1]);
                        Seed_point2=[Seed_point2;Weight_Pnt];
                    else
                        y1=round((D*d(i,1)-x(i)*y1my2-C)/x2mx1);
                        yy2=[x(i),y1;yy2];
                        Seed_point2=[Seed_point2;[x(i),y1]];
                        
                    end
                end
                
                if  (m<=5)  &   (m>-4)
                    Seed_point=intersect(Seed_point2,Seed_point1,'rows');
                    yy=intersect(yy1,yy2,'rows');
                else
                    Seed_point=Seed_point2;
                    yy=yy2;
                end
                % find NaN and Inf
                 X_temp=Seed_point(:,1);
                idx_inf=find(X_temp==Inf);
                Seed_point(idx_inf,:)=[];
                 X_temp=Seed_point(:,1);
                idx_nan=find(isnan(X_temp)==1);
                Seed_point(idx_nan,:)=[];
                Y_temp=Seed_point(:,2);
                idx2_inf=find(Y_temp==Inf);
                Seed_point(idx2_inf,:)=[];
                 Y_temp=Seed_point(:,2);
                idx2_nan=find(isnan(Y_temp)==1);
                Seed_point(idx2_nan,:)=[];
                 X_temp=Seed_point(:,1);
                 Y_temp=Seed_point(:,2);
                

                idx=randperm(length(Seed_point));
               
                if length(X_temp)>20
                    N_Seed=length(X_temp);
                else
                    N_Seed=length(X_temp);
                end
                Seed1=X_temp(idx(1:N_Seed));
                Seed2=Y_temp(idx(1:N_Seed));
                if isempty(Seed1)| isempty(Seed2)
                   CurveEdge=CurveEdge-1;
                else
                 Point_Between{CurveEdge}=yy;
                 Seed_point_6{CurveEdge}=[Seed1 Seed2];
                end
              
        end
    end
end


toc
Segment.seedpoint=Seed_point_6;

Segment.hypo=recog_result.hypo_mask;
Segment.edge=Prob_Carried_Object;
Segment.PsnVeiw=recog_result.PsnVeiw;

%save(['/home/farnoosh/Projects/CO_Detector/Hypo/Hypo_iLid/',Img,'_',num2str(PN), '.mat'],'Segment');

%drawedgelist(Seed_point_6,Point_Between,New_edgelist,edgelist,Curve_List, size(Obj_Msk), 1, 'rand', 2); 
Hypo_Prob=Segment;
% 
