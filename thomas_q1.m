% Defining the input matrices
a=[0;2;3;6;1];
b=[4;4;-9;8;7];
c=[-2;1;2;1;0];
d=[1;2;-3;4;5];
n=5;
x=thomas_q1(a,b,c,d,n);

disp("a=")
disp(a);
disp("b=")
disp(b);
disp("c=")
disp(c);
disp("d=")
disp(d);
disp("The solution to the tri-diagonal system");
disp(x);



function x = thomas_q1(a,b,c,d,n)

% Matrix a contains terms below the main diagonal 
% Matrix b contains terms on the main diagonal 
% Matrix c contains terms above the main diagonal 
% Matrix d is the RHS of the system

%% Reduce the tri-diagonal matrix by row transformations

for i=1:n-1
  
 % Normalize by dividing with diagonal element   
    c(i,1)=c(i,1)/b(i,1);
    d(i,1)=d(i,1)/b(i,1);
    b(i,1)=1;
 % Convert the matrix A to the Echelon form    
    alpha=a(i+1,1);
    a(i+1,1)=0;
    b(i+1,1)=b(i+1,1)-alpha*c(i,1);
    d(i+1,1)=d(i+1,1)-alpha*d(i,1);
end

%% Initialize the vector of solutions 
x=zeros(n,1); 

%% Back-substitution to find out the solution
x(n,1)=d(n,1)/b(n,1);

for i=n-1:-1:1
    x(i,1)=d(i,1)-x(i+1,1)*c(i,1);
end

end

