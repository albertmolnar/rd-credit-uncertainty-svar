function c = setf(c, varargin)
% Returns cfg c with field/value pairs overwritten: setf(c, 'p', 2, 'Hrestr', 0).
for i = 1:2:numel(varargin)
    c.(varargin{i}) = varargin{i+1};
end
end
