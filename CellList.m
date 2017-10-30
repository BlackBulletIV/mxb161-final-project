classdef CellList < handle
    properties
        First = false
        Last = false
        Length = 0
    end
    
    methods
        function cell = add(self, x, y)
            cell = Cell(x, y, self.Last, false, self);

            if islogical(self.First)
                self.First = cell;
            end
            
            if ~islogical(self.Last)
                self.Last.Next = cell;
            end
            
            self.Last = cell;
            self.Length = self.Length + 1;
        end
        
        function remove(self, cell)
            if ~islogical(cell.List)
                if ~islogical(cell.Prev)
                    if islogical(cell.Next)
                        self.Last = cell.Prev;
                        cell.Prev.Next = false;
                    else
                        cell.Prev.Next = cell.Next;
                        cell.Next.Prev = cell.Prev;
                    end
                elseif ~islogical(cell.Next)
                    self.First = cell.Next;
                    cell.Next.Prev = false;
                else
                    self.First = false;
                    self.Last = false;
                end

                self.Length = self.Length - 1;
                cell.List = false;
            end
        end            
    end 
end
