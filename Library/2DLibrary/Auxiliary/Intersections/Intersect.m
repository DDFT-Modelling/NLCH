function data = Intersect(MainShape,SecondShape)

    intersectFnStr = ['Intersect_' class(SecondShape)];

    if(exist(intersectFnStr,'file')==2) % is a Matlab file
        intersectFn = str2func(intersectFnStr);
        data = intersectFn(MainShape,SecondShape);
    else
        exc = MException('Intersect','case not implemented');
        throw(exc);
    end

end