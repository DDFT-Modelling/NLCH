classdef ArbitraryIntersection < handle

    properties
        Pts
        Int
        Area
        Origin = [0;0];
    end
   
    methods
        function this = ArbitraryIntersection(Shapes)
            if(nargin == 0)
                nShapes = 0;
            else
                nShapes = length(Shapes);
            end
            
            y1_kv = [];
            y2_kv = [];
            this.Int = [];
            this.Area = 0;
            
            for iShape = 1:nShapes
                shapePts = Shapes(iShape).Shape.GetCartPts;
                y1_kv = [y1_kv ; shapePts.y1_kv];
                y2_kv = [y2_kv ; shapePts.y2_kv];
                [int,area] = Shapes(iShape).Shape.ComputeIntegrationVector;                    
                this.Int = [this.Int , int];
                this.Area = this.Area + area;
            end
            
            this.Pts.y1_kv = y1_kv;
            this.Pts.y2_kv = y2_kv;
            
            if(nargin>0)
                this.Origin = Shapes(1).Shape.Origin;
            end
        end
        
        
        function ptsCart = GetCartPts(this)            
            ptsCart       = this.Pts;
        end                 
        
        
        function [int,area] = ComputeIntegrationVector(this)                       
            int = this.Int;
            area = this.Area;
        end
        
    end
    
    
    
end