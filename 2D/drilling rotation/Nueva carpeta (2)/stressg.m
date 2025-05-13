function [stressg] = stressg(forcm,ro)

% ********   forces internes dans le repere global   ********

      stressg(1,1)=forcm(1,1)*ro(1,1)^2+forcm(1,2)*ro(2,1)^2+2*forcm(1,3)*ro(1,1)*ro(2,1);
      stressg(1,2)=forcm(1,1)*ro(1,2)^2+forcm(1,2)*ro(2,2)^2+2*forcm(1,3)*ro(1,2)*ro(2,2);
      stressg(1,3)=forcm(1,2)*ro(1,1)*ro(2,1)-forcm(1,1)*ro(1,1)*ro(2,1)+forcm(1,3)*(ro(1,1)^2-ro(2,1)^2);
      stressg(2,1)=forcm(2,1)*ro(1,1)^2+forcm(2,2)*ro(2,1)^2+2*forcm(2,3)*ro(1,1)*ro(2,1);
      stressg(2,2)=forcm(2,1)*ro(1,2)^2+forcm(2,2)*ro(2,2)^2+2*forcm(2,3)*ro(1,2)*ro(2,2);
      stressg(2,3)=forcm(2,2)*ro(1,1)*ro(2,1)-forcm(2,1)*ro(1,1)*ro(2,1)+forcm(2,3)*(ro(1,1)^2-ro(2,1)^2);
      stressg(3,1)=forcm(3,1)*ro(1,1)^2+forcm(3,2)*ro(2,1)^2+2*forcm(3,3)*ro(1,1)*ro(2,1);
      stressg(3,2)=forcm(3,1)*ro(1,2)^2+forcm(3,2)*ro(2,2)^2+2*forcm(3,3)*ro(1,2)*ro(2,2);
      stressg(3,3)=forcm(3,2)*ro(1,1)*ro(2,1)-forcm(3,1)*ro(1,1)*ro(2,1)+forcm(3,3)*(ro(1,1)^2-ro(2,1)^2);

end