% This script will take input from the user about the beam and 
% shows its bending moment and shear force diagram
% This script is only for simply supported beam

%input part
clc; clearvars;


disp('We are dealing with a simply supported beam with supports at its ends')
disp(' ')
len = input('Enter the length of beam in meters : ');
disp(' ');

disp('Enter the loads from left to right ');
disp(' ');

%point loads
num = input('Enter the number of point loads : ');
disp(' ');

if num > 0
    point_loads = zeros(1,num);
    pos_point = zeros(1,num);


    for i = 1:num
    
        x = ['Enter the values of load ', num2str(i), ' in kN : '];
        point_loads(1,i) = input(x);
        
    
        x = ['Enter the position of load ', num2str(i), ' (from left end) in meters : '];
        pos_point(1,i) = input(x);
        disp(' ');
    
        while pos_point(1,i) > len || pos_point(1,i) < 0 
            disp('Invalid load position');
            fprintf('Enter the value between 0 and %d \n',len);
            disp(' ');
            pos_point(1,i) = input('Enter the position of load (from left end) in meters : ');
            disp(' ');
        end
    end
end

% uniformly distributed load
u_num = input('Enter number of uniformly distributed loads : ');
disp(' ');

if u_num > 0
    
    u_load = zeros(1,u_num);
    u_pos = zeros(2,u_num);
    % first row for starting point second row for ending point
    
    for i = 1:u_num
        
        u_load(1, i) = input('Enter the value of load per unit length in kN/m : ');
        u_pos(1, i) = input('Enter the starting point of force from left end in m : ');
        % retake value if invlid value is entered
        while u_pos(1, i) < 0 || u_pos(1, i) > len
            disp(' ');
            str = ['Enter a valid value between 0 and ' num2str(len) ' '];
            disp(str);
            disp(' ');
            u_pos(1, i) = input('Enter the starting point of force from left end in m : ');
        end
            
        u_pos(2, i) = input('Enter the ending point of force from left end in m : ');
        % retake value if invalid input is entered
        while u_pos(2, i) < u_pos(1, i) || u_pos(2, i) > len
            disp(' ');
            str = ['Enter valid value between ' num2str(u_pos(1,i)) ' and ' num2str(len) ' '];
            disp(str);
            disp(' ');
            u_pos(2, i) = input('Enter the ending point of force from left end in m : ');
        end
        disp(' ');
    end
    
end

% additional moments

m_num = input('Enter the number of moments on beam : ');
disp(' ');

m_val = zeros(1, m_num);
m_pos = zeros(1, m_num);

for i = 1:m_num
    str = ['Enter value of moment ' num2str(i) ' : '];
    m_val(1, i) = input(str);
    str = ['Enter the position of moment ' num2str(i) ' from left end : '];
    m_pos(1, i) = input(str);
    disp(' ');
    
    while m_pos(1, i) < 0 || m_pos(1, i) > len
        str = ['Enter a valid value between 0 and ' len];
        disp(str);
        disp(' ');
        
        str = ['Enter the position of moment ' num2str(i) ' from left end : '];
        m_pos(1, i) = input(str);
        disp(' ');
    end
end
        
        
% calculating part
x = 0:1/1000:len;

% calculate support forces
a_y = 0;
b_y = 0;

% point loads
for i = 1:num
    
    a_y = a_y + point_loads(1,i)*(len - pos_point(1,i));
    b_y = b_y + point_loads(1,i)*(pos_point(1,i));
    
end

% uniformly distributed forces
for i = 1:u_num
    
    a_y = a_y + u_load(1, i)*(u_pos(2, i) - u_pos(1, i))*(len - (u_pos(2, i) + u_pos(1, i))/2);
    b_y = b_y + u_load(1, i)*(u_pos(2, i) - u_pos(1, i))*((u_pos(2, i) + u_pos(1, i))/2);
    
end

% additional moments
for i = 1:m_num
    
    a_y = a_y + m_val(1, i);
    b_y = b_y + m_val(1, i);
    
end

a_y = a_y/len ;
b_y = b_y/len ;

str = ['Reaction forces are :' newline 'A_y = ' num2str(a_y) 'kN' newline 'B_y = ' num2str(b_y) 'kN'];
disp(str);
%fprintf("Reaction forces are : \n A_y = %.2d kN \n B_y = %.2d kN \n", a_y, b_y);

% calculate 'V' shear force and 'M'

V = -1*a_y + 0*x;
M = a_y*x;

% point loads
for i = 1:num
    % Using Macaulay's method for point loads
    V = V + (x > pos_point(1,i))*point_loads(1,i);
    M = M - (x > pos_point(1,i)).*(point_loads(1,i).*(x - pos_point(1,i)));
end

% u_loads
for i = 1:u_num
     
    V = V + ((x > u_pos(1, i)) & (x < u_pos(2, i))).*u_load(i).*(x - u_pos(1, i)) + (x >= u_pos(2, i))*u_load(i)*(u_pos(2,i) - u_pos(1,i));
    M = M - ((x > u_pos(1, i)) & (x < u_pos(2, i))).*u_load(i).*((x - u_pos(1, i)).^(2) /2) - (x >= u_pos(2, i)).*u_load(i).*(u_pos(2,i) - u_pos(1,i)).*(x - (u_pos(2, i) + u_pos(1, i))/2); 

end

% additional moments
for i = 1:m_num
    M = M - (x >= m_pos(1, i))*m_val(1, i);
end

%plotting diagram
ax1 = subplot(3,1,1);
hold on
line([0 len], [0 0],'Color' , 'k', 'LineWidth', 2);
axis([-1 len+1 -2 2]);

% arrows for a_y and b_y
quiver(0, -1, 0, 1, 1, 'Color', 'b', 'Linewidth', 1);     % a_y
str = [num2str(a_y) ' kN'];
text(0, -1.2, str);

str = [num2str(b_y) ' kN'];
text(len, -1.2, str); 
quiver(len, -1, 0, 1, 1, 'Color', 'b', 'Linewidth', 1);     % b_y

% arrows for point loads
for i = 1:num
    if point_loads(1, i) >= 0 
        quiver(pos_point(1, i), 1.3, 0, -1.3, 1, 'Color', 'r', 'Linewidth', 1);     % down arrow
        str = [num2str(point_loads(1, i)) ' kN'];
        text(pos_point(1, i), 1.5, str);
        
    else
        quiver(pos_point(1, i), -1.3, 0, 1.3, 1, 'Color', 'r', 'Linewidth', 1);    % up arrow
        str = [num2str(point_loads(1, i)) ' kN'];
        text(pos_point(1, i), -1.5, str);
    end
end

% arrows for uniformly distributed loads
for i = 1:u_num
    for j = u_pos(1, i):0.25:u_pos(2, i)
        if u_load(1, i) >= 0
            quiver(j, 1, 0, -1, 1, 'Color', 'b', 'Linewidth', 1);     % down arrows
        else
            quiver(j,-1, 0, 1, 1, 'Color', 'b', 'Linewidth', 1);       % up arrows
        end
    end
    
    if u_load(1, i) >= 0   % joining line and text
        line([u_pos(1, i) u_pos(2, i)], [1 1], 'Color', 'b');
        str = [num2str(u_load(1,i)) ' kN/m'];
        text((u_pos(1,i) + u_pos(2,i))/2, 1.2, str);
    else
        line([u_pos(1, i) u_pos(2, i)], [-1 -1], 'Color', 'b');
        str = [num2str(u_load(1,i)) 'kN/m'];
        text((u_pos(1,i) + u_pos(2,i))/2, -1.2, str);
    end
end

% stars for additional moments

for i = 1:m_num
    plot(ax1, m_pos(1, i), 0, '*g');
    str = [num2str(m_val(1, i)) 'kN-m'];
    text(m_pos(1, i), 0.5, str);
end

%             down arrow  //////////////
%quiver([x], [y], [0], [-1], 0.5);

%             up arrow    ///////////////
%quiver([x], [y], [0], [1], 0.5);



hold(ax1, 'off');

%plotting shear force diagram
ax2 = subplot(3,1,2);
plot(ax2,x,V);
title('Shear force');
hold(ax2, 'on');
grid on

v_points = zeros(1, num + 2*u_num + m_num);
m_points = zeros(1, num + 2*u_num + m_num);
pos = zeros(1, num + 2*u_num + m_num);

for i = 1:(num + 2*u_num + m_num)
    if i <= num
        v_points(1, i) = V(1, 1+pos_point(1, i)*1000);
        m_points(1, i) = M(1, 1+pos_point(1, i)*1000);
        pos(1, i) = pos_point(1, i);
    elseif i<= num + u_num
        v_points(1 ,i) = V(1, 1+u_pos(1, i - num)*1000);
        m_points(1, i) = M(1, 1+u_pos(1, i - num)*1000);
        pos(1, i) = u_pos(1, i - num);
    elseif i<= num + 2*u_num
        v_points(1, i) = V(1, 1+u_pos(2, i - num - u_num)*1000);
        m_points(1, i) = M(1, 1+u_pos(2, i - num - u_num)*1000);
        pos(1, i) = u_pos(2, i - num - u_num);
    else 
        v_points(1, i) = V(1, 1+m_pos(1, i - num - 2*u_num)*1000);
        m_points(1, i) = M(1, 1+m_pos(1, i - num - 2*u_num)*1000);
        pos(1, i) = m_pos(1, i - num - 2*u_num);
    end
        
end


plot(ax2, pos, v_points,'or','MarkerFaceColor','y');
line([0 len], [0 0], 'Color', 'k','LineWidth',1);  %axis
line([0, 0], [0, V(1,1)]); % starting line of shear force
line([len, len], [0, V(1,len*1000)]); % ending line of shear force
% fill(x, V, 'b');
v_max = max(V);
v_min = min(V);
axis([-1 len+1 v_min*(13/10) v_max*(13/10)]);
hold(ax2, 'off');


% plotting bending moment diagram
ax3 = subplot(3,1,3);
plot(ax3,x,M);
grid on
title('Bending Moment');
hold(ax3, 'on');
line([0 len], [0 0], 'Color', 'k','LineWidth',1);
% fill(x, M, 'b');
plot(ax3, pos, m_points, 'or', 'MarkerFaceColor','y');
line([0 len], [0 0], 'Color', 'k','LineWidth',1); % axis
m_max = max(M);
m_min = min(M);
axis([-1 len+1 m_min*(13/10) m_max*(13/10)]);
hold(ax3, 'off');



    




