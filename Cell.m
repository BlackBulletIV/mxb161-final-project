classdef Cell < handle
    properties
        X
        Y
        DepartX
        DepartY
        Course
        CanLand = 0
        Prev = false
        Next = false
        List = false
        PX = 0
        PY = 0
        Times = 0
    end
    
    methods
        function cell = Cell(x, y, prev, next, list)
            cell.X = x;
            cell.Y = y;
            cell.Prev = prev;
            cell.Next = next;
            cell.List = list;
        end
    end
end
