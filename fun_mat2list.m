function list_val = fun_mat2list(mat)

l           = size(mat);

a           = string(num2str(mat,'%+32.32f\\'));
ind         = strfind(a,'\');
b           = char(a);
ind         = ind{1};
b(:,ind)    = ',';
lef_bracket = repmat('[',l(1),1);
rgt_bracket = repmat(']',l(1),1);
cma_val     = repmat(',',l(1),1);
b           = [lef_bracket,b,rgt_bracket,cma_val];
list_val    = ['[',reshape(b',1,size(b,1) * size(b,2)),']'];


end