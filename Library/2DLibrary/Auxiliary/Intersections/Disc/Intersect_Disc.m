function dataDisc = Intersect_Disc(MainShape,discShape)

    intersectFnStr = ['Intersect_Disc_' class(MainShape)];

    if(exist(intersectFnStr,'file')==2) % is a Matlab file
        intersectFn = str2func(intersectFnStr);
        area = intersectFn(discShape,MainShape);
    else
        exc = MException('IntersectDisc',['case' intersectFnStr 'not implemented']);
        throw(exc);
    end

    y10 = discShape.Origin(1);   
	y20 = discShape.Origin(2);       

    % get polar coordinates for, e.g., FMT kernels
    dataDisc.pts       = area.GetCartPts();
    ptsLoc.y1_kv       = dataDisc.pts.y1_kv - y10;
    ptsLoc.y2_kv       = dataDisc.pts.y2_kv - y20;
    dataDisc.ptsPolLoc = Cart2PolPts(ptsLoc);

    % deal with intersections at infinity
    if(y10 == inf)
        dataDisc.pts.y1_kv = dataDisc.pts.y1_kv + y10;
    end
    if(y20 == inf)
        dataDisc.pts.y2_kv = dataDisc.pts.y2_kv + y20;
    end
    
    [dataDisc.int,dataDisc.area]     = area.ComputeIntegrationVector();                                   
    
end 

%     if(isa(MainShape,'HalfSpace'))   
%         area = Intersect_Disc_HalfSpace(discShape,MainShape);
%         
%     elseif(isa(MainShape,'Box'))
%         area           = Intersect_Disc_Box(discShape,MainShape);                
% 
%     elseif(isa(MainShape,'PeriodicBox'))
%         area           = Intersect_Disc_PeriodicBox(discShape,MainShape);                
%           
%     elseif(isa(MainShape,'InfCapillary') && (MainShape.y2Max - MainShape.y2Min)>= 2*discShape.R)   
% 
%         area = Intersect_Disc_InfCapillary(discShape,MainShape);
%     
%     else
%         exc = MException('Intersect_Disc','case not implemented');
%         throw(exc);                
%     end        
