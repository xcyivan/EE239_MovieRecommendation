Data_raw = importdata('u.data');
userid = Data_raw(:,1);
itemid = Data_raw(:,2);
rating = Data_raw(:,3);
row = max(userid);
col = max(itemid);
R = zeros(row,col);
W = zeros(row,col);
 for i = 1:size(userid)
     R (userid(i),itemid(i)) = rating(i); 
 end
 W((R>0))=1;
%/////////Problem 1
%dicompose R into U,V, k =10
%Error matrix, error1 = 5.4620e+04
 [U1,V1] = wnmfrule (R,10);
 E1 = W.*((R-U1*V1).^2);
 error1 = sum(E1(:));
 
 %repeat for k = 50, error2 = 1.9493e+04
 [U2,V2] = wnmfrule (R, 50);
 E2 = W.*((R-U2*V2).^2);
 error2 = sum(E2(:));
 
 %repeat for k = 100, error3 = 6.2705e+03 
 [U3,V3] = wnmfrule (R, 100);
 E3 = W.*((R-U3*V3).^2);
 error3 = sum(E3(:));
 

 %/////////Problem 2
 %////System2
 %dicompose R into U,V, k =10
 %Error matrix, error4 = 1.6260e-06
 [U4,V4] = wnmfrule2 (R,10);
 E4 = R.*((W-U4*V4).^2);
 error4 = sum(E4(:));
 
 %repeat for k = 50, error5 = 6.0423e-05
 [U5,V5] = wnmfrule2 (R, 50);
 E5 = R.*((W-U5*V5).^2);
 error5 = sum(E5(:));
 
 %repeat for k = 100, error6 = 1.1560e-04 
 [U6,V6] = wnmfrule2 (R, 100);
 E6 = R.*((W-U6*V6).^2);
 error6 = sum(E6(:));
 
 %///////Problem 3
known_indices = find(R~=0);
N = length(known_indices);
random_indices = randperm(N);


%test1
test1_indices = known_indices(random_indices(1+(1-1)*floor(N/10):1*floor(N/10)));
R1 = R;
R1(test1_indices) = 0;
[A1,Y1] = wnmfrule(R1,100);
R1_rec = A1*Y1;
test1_rating = R(test1_indices);
pred1_rating = R2_rec(test1_indices);
errdis1=sum(abs(pred1_rating-test1_rating));
precision1 = length(find(test1_rating>3 & pred1_rating>3))/length(find(pred1_rating>3));
recall1 = length(find(test1_rating>3 & pred1_rating>3)) / length(find(test1_rating>3));


%test2
test2_indices = known_indices(random_indices(1+(2-1)*floor(N/10):2*floor(N/10)));
R2 = R;
R2(test2_indices) = 0;
[A2,Y2] = wnmfrule(R2,100);
R2_rec = A2*Y2;
test2_rating = R(test2_indices);
pred2_rating = R2_rec(test2_indices);
errdis2=sum(abs(pred2_rating-test2_rating));
precision2 =  length(find(test2_rating>3 & pred2_rating>3))/length(find(pred2_rating>3));
recall2 = length(find(test2_rating>3 & pred2_rating>3)) / length(find(test2_rating>3));

%test3
test3_indices = known_indices(random_indices(1+(3-1)*floor(N/10):3*floor(N/10)));
R3 = R;
R3(test3_indices) = 0;
[A3,Y3] = wnmfrule(R3,100);
R3_rec = A3*Y3;
test3_rating = R(test3_indices);
pred3_rating = R3_rec(test3_indices);
errdis3=sum(abs(pred3_rating-test3_rating));
precision3 =  length(find(test3_rating>3 & pred3_rating>3))/length(find(pred3_rating>3));
recall3 = length(find(test3_rating>3 & pred3_rating>3)) / length(find(test3_rating>3));



%test4
test4_indices = known_indices(random_indices(1+(4-1)*floor(N/10):4*floor(N/10)));
R4 = R;
R4(test4_indices) = 0;
[A4,Y4] = wnmfrule(R4,100);
R4_rec = A4*Y4;
test4_rating = R(test4_indices);
pred4_rating = R4_rec(test4_indices);
errdis4=sum(abs(pred4_rating-test4_rating));
precision4 =  length(find(test4_rating>3 & pred4_rating>3))/length(find(pred4_rating>3));
recall4 = length(find(test4_rating>3 & pred4_rating>3)) / length(find(test4_rating>3));



%test5
test5_indices = known_indices(random_indices(1+(5-1)*floor(N/10):5*floor(N/10)));
R5 = R;
R5(test5_indices) = 0;
[A5,Y5] = wnmfrule(R5,100);
R5_rec = A5*Y5;
test5_rating = R(test5_indices);
pred5_rating = R5_rec(test5_indices);
errdis5=sum(abs(pred5_rating-test5_rating));
precision5 =  length(find(test5_rating>3 & pred5_rating>3))/length(find(pred5_rating>3));
recall5 = length(find(test5_rating>3 & pred5_rating>3)) / length(find(test5_rating>3));



%test6
test6_indices = known_indices(random_indices(1+(6-1)*floor(N/10):6*floor(N/10)));
R6 = R;
R6(test6_indices) = 0;
[A6,Y6] = wnmfrule(R6,100);
R6_rec = A6*Y6;
test6_rating = R(test6_indices);
pred6_rating = R6_rec(test6_indices);
errdis6=sum(abs(pred6_rating-test6_rating));
precision6 =  length(find(test6_rating>3 & pred6_rating>3))/length(find(pred6_rating>3));
recall6 = length(find(test6_rating>3 & pred6_rating>3)) / length(find(test6_rating>3));



%test7
test7_indices = known_indices(random_indices(1+(7-1)*floor(N/10):7*floor(N/10)));
R7 = R;
R7(test7_indices) = 0;
[A7,Y7] = wnmfrule(R7,100);
R7_rec = A7*Y7;
test7_rating = R(test7_indices);
pred7_rating = R7_rec(test7_indices);
errdis7=sum(abs(pred7_rating-test7_rating));
precision7 =  length(find(test7_rating>3 & pred7_rating>3))/length(find(pred7_rating>3));
recall7 = length(find(test7_rating>3 & pred7_rating>3)) / length(find(test7_rating>3));



%test8
test8_indices = known_indices(random_indices(1+(8-1)*floor(N/10):8*floor(N/10)));
R8 = R;
R8(test8_indices) = 0;
[A8,Y8] = wnmfrule(R8,100);
R8_rec = A8*Y8;
test8_rating = R(test8_indices);
pred8_rating = R8_rec(test8_indices);
errdis8=sum(abs(pred8_rating-test8_rating));
precision8 =  length(find(test8_rating>3 & pred8_rating>3))/length(find(pred8_rating>3));
recall8 = length(find(test8_rating>3 & pred8_rating>3)) / length(find(test8_rating>3));



%test9
test9_indices = known_indices(random_indices(1+(9-1)*floor(N/10):9*floor(N/10)));
R9 = R;
R9(test9_indices) = 0;
[A9,Y9] = wnmfrule(R9,100);
R9_rec = A9*Y9;
test9_rating = R(test9_indices);
pred9_rating = R9_rec(test9_indices);
errdis9=sum(abs(pred9_rating-test9_rating));
precision9 =  length(find(test9_rating>3 & pred9_rating>3))/length(find(pred9_rating>3));
recall9 = length(find(test9_rating>3 & pred9_rating>3)) / length(find(test9_rating>3));



%test10
test10_indices = known_indices(random_indices(1+(10-1)*floor(N/10):10*floor(N/10)));
R10 = R;
R10(test10_indices) = 0;
[A10,Y10] = wnmfrule(R10,100);
R10_rec = A10*Y10;
test10_rating = R(test10_indices);
pred10_rating = R10_rec(test10_indices);
errdis10=sum(abs(pred10_rating-test10_rating));
precision10 =  length(find(test10_rating>3 & pred10_rating>3))/length(find(pred10_rating>3));
recall10 = length(find(test10_rating>3 & pred10_rating>3)) / length(find(test10_rating>3));



%System 2
%test11
[A11,Y11] = wnmfrule2(R1,10);
R11_rec = (A11*Y11).*R;
test11_rating = test1_rating;
pred11_rating = R11_rec(test1_indices);
errdis11 = sum(abs(pred11_rating-test11_rating));
% precision11 =  length(find(test11_rating>3 & pred11_rating>3))/length(find(pred11_rating>3));
% recall11 = length(find(test11_rating>3 & pred11_rating>3)) / length(find(test11_rating>3));



%test12
[A12,Y12] = wnmfrule2(R2,10);
R12_rec = (A12*Y12).*R;
test12_rating = test2_rating;
pred12_rating = R12_rec(test2_indices);
errdis12 = sum(abs(pred12_rating-test12_rating));
% precision12 =  length(find(test12_rating>3 & pred12_rating>3))/length(find(pred12_rating>3));
% recall12 = length(find(test12_rating>3 & pred12_rating>3)) / length(find(test12_rating>3));


%test13
[A13,Y13] = wnmfrule2(R3,10);
R13_rec = (A13*Y13).*R;
test13_rating = test3_rating;
pred13_rating = R13_rec(test3_indices);
errdis13 = sum(abs(pred13_rating-test13_rating));
% precision13 =  length(find(test13_rating>3 & pred13_rating>3))/length(find(pred13_rating>3));
% recall13 = length(find(test13_rating>3 & pred13_rating>3)) / length(find(test13_rating>3));



%test14
[A14,Y14] = wnmfrule2(R4,10);
R14_rec = (A14*Y14).*R;
test14_rating = test4_rating;
pred14_rating = R14_rec(test4_indices);
errdis14 = sum(abs(pred14_rating-test14_rating));
% precision14 =  length(find(test14_rating>3 & pred14_rating>3))/length(find(pred14_rating>3));
% recall14 = length(find(test14_rating>3 & pred14_rating>3)) / length(find(test14_rating>3));



%test15
[A15,Y15] = wnmfrule2(R5,10);
R15_rec = (A15*Y15).*R;
test15_rating = test5_rating;
pred15_rating = R15_rec(test5_indices);
errdis15 = sum(abs(pred15_rating-test15_rating));
% precision15 =  length(find(test15_rating>3 & pred15_rating>3))/length(find(pred15_rating>3));
% recall15 = length(find(test15_rating>3 & pred15_rating>3)) / length(find(test15_rating>3));



%test16
[A16,Y16] = wnmfrule2(R6,10);
R16_rec = (A16*Y16).*R;
test16_rating = test6_rating;
pred16_rating = R16_rec(test6_indices);
errdis16 = sum(abs(pred16_rating-test16_rating));
% precision16 =  length(find(test16_rating>3 & pred16_rating>3))/length(find(pred16_rating>3));
% recall16 = length(find(test16_rating>3 & pred16_rating>3)) / length(find(test16_rating>3));



%test17
[A17,Y17] = wnmfrule2(R7,10);
R17_rec = (A17*Y17).*R;
test17_rating = test7_rating;
pred17_rating = R17_rec(test7_indices);
errdis17 = sum(abs(pred17_rating-test17_rating));
% precision17 =  length(find(test17_rating>3 & pred17_rating>3))/length(find(pred17_rating>3));
% recall17 = length(find(test17_rating>3 & pred17_rating>3)) / length(find(test17_rating>3));



%test18
[A18,Y18] = wnmfrule2(R8,10);
R18_rec = (A18*Y18).*R;
test18_rating = test8_rating;
pred18_rating = R18_rec(test8_indices);
errdis18 = sum(abs(pred18_rating-test18_rating));
% precision18 =  length(find(test18_rating>3 & pred18_rating>3))/length(find(pred18_rating>3));
% recall18 = length(find(test18_rating>3 & pred18_rating>3)) / length(find(test18_rating>3));



%test19
[A19,Y19] = wnmfrule2(R9,10);
R19_rec = (A19*Y19).*R;
test19_rating = test9_rating;
pred19_rating = R19_rec(test9_indices);
errdis19 = sum(abs(pred19_rating-test19_rating));
% precision19 =  length(find(test19_rating>3 & pred19_rating>3))/length(find(pred19_rating>3));
% recall19 = length(find(test19_rating>3 & pred19_rating>3)) / length(find(test19_rating>3));



%test20
[A20,Y20] = wnmfrule2(R10,10);
R20_rec = (A20*Y20).*R;
test20_rating = test10_rating;
pred20_rating = R20_rec(test10_indices);
errdis20 = sum(abs(pred20_rating-test20_rating));
% precision20 =  length(find(test20_rating>3 & pred20_rating>3))/length(find(pred20_rating>3));
% recall20 = length(find(test20_rating>3 & pred20_rating>3)) / length(find(test20_rating>3));



syserr1=(errdis1+errdis2+errdis3+errdis4+errdis5+errdis6+errdis7+errdis8+errdis9+errdis10)/10;
syserr2=(errdis11+errdis12+errdis13+errdis14+errdis15+errdis16+errdis17+errdis18+errdis19+errdis20)/10;

%Problem 4
%choose test 7 and test 17 because they give least error distance
precision_sys1=zeros(51);
recall_sys1=zeros(51);
precision_sys2=zeros(51);
recall_sys2=zeros(51);
temp17 = A17*Y17;
temp17_2 = A17_2*Y17_2;
pred17_raw = temp17(test7_indices);
pred17_2_raw= temp17_2(test7_indices);
j=1;
for i = 0:0.1:5
    precision_sys1(j)=length(find(test7_rating>3 & pred7_rating>i))/length(find(pred7_rating>i));
    recall_sys1(j) = length(find(test7_rating>3 & pred7_rating>i)) / length(find(test7_rating>3));
    precision_sys2(j) =length(find(test17_rating>3 & pred17_rating>i))/length(find(pred17_rating>i));
    recall_sys2(j) = length(find(test17_rating>3 & pred17_rating>i)) / length(find(test17_rating>3));
    j=j+1;
end
j=1;
for i=0:0.02:1
    precision_sys3(j) =length(find(test17_rating>3 & pred17_raw>i))/length(find(pred17_raw>i));
    precision_sys3_2(j)=length(find(test17_rating>3 & pred17_2_raw>i))/length(find(pred17_2_raw>i));
    recall_sys3(j) = length(find(test17_rating>3 & pred17_raw>i)) / length(find(test17_rating>3));
    recall_sys3_2(j)=length(find(test17_rating>3 & pred17_2_raw>i))/length(find(test17_rating>3));
    j=j+1;
end
x=(0:0.1:5);
xx=(0:0.02:1);
subplot(2,1,1);
plot(x,precision_sys1);
xlabel('threshold');
ylabel('precision');
title('Sytem 1 : precision');
axis([0,5,0.5,0.75]);
subplot(2,1,2);
plot(x,recall_sys1);
xlabel('threshold');
ylabel('recall');
title('Sytem 1 : recall');
axis([0,5,0,1]);

%R_rec .R, not actually meet the requirement but make more sense
figure;
subplot(2,1,1);
plot(x,precision_sys2);
xlabel('threshold');
ylabel('precision');
title('Sytem 2 .*R: precision');
subplot(2,1,2);
plot(x,recall_sys2);
xlabel('threshold');
ylabel('recall');
title('System 2 .*R: recall');

%actually the figure of system 2 that meets the requirement
figure;
subplot(2,1,1);
plot(xx,precision_sys3);
xlabel('threshold');
ylabel('precision');
title('System 2 : precision');
subplot(2,1,2);
plot(xx,recall_sys3);
xlabel('threshold');
ylabel('recall');
title('System 2 : recall');

%using method in Problem 5
figure;
subplot(2,1,1);
plot(xx,precision_sys3_2);
xlabel('threshold');
ylabel('precision');
title('System 2 regularized: precision');
subplot(2,1,2);
plot(xx,recall_sys3_2);
xlabel('threshold');
ylabel('recall');
title('System 2 regularized: recall');
%////Problem 5, 
%Compared with System 1, take k = 100
%lamda = 0.01
 [U7,V7] = wnmfrule5 (0.01,R,100);
 E7 = W.*((R-U7*V7).^2);
 error7 = sum(E7(:));
 
 %lamda = 0.1
 [U8,V8] = wnmfrule5 (0.1,R,100);
 E8 = W.*((R-U8*V8).^2);
 error8 = sum(E8(:));
 
 %lamda = 1
 [U9,V9] = wnmfrule5 (1,R,100);
 E9 = W.*((R-U9*V9).^2);
 error9 = sum(E9(:));
 
 %Compared with System 2, take k =10;
 %lamda = 0.01
 [U10,V10] = wnmfrule6 (0.01,R,10);
 E10 = R.*((W-U10*V10).^2);
 error10 = sum(E10(:));
 
  %lamda = 0.1
 [U11,V11] = wnmfrule6 (0.1,R,10);
 E11 = R.*((W-U11*V11).^2);
 error11 = sum(E11(:));
 
  %lamda = 1
 [U12,V12] = wnmfrule6 (1,R,10);
 E12 = R.*((W-U12*V12).^2);
 error12 = sum(E12(:));

 
 %update problem 5 System2, take lamda = 0.01 cz it yieds the least error
[A17_2,Y17_2] = wnmfrule6(0.01,R7,10);
R17_2_rec = (A17_2*Y17_2).*R;
test17_2_rating = test7_rating;
pred17_2_rating = R17_2_rec(test7_indices);
 

%Problem 6
%Problem 6 update
R_evaluate_0 = zeros(row,col); % those predict using 0-1 matrix, as the requirement says
R_evaluate = zeros(row,col); % new algorhtms with A*Y .* Rating // although not fit the requirement but give great performance
R_evaluate_2 = zeros(row,col); % using method in Problem 5, performance worse than .*Rating algorithm
R_origin = zeros(row,col);
R_rec_0 = A17*Y17;
R_rec = A17*Y17.*R;
R_rec_2 = A17_2 * Y17_2;
R_evaluete_0(test7_indices) = R_rec_0(test7_indices);
R_evaluate(test7_indices) = R_rec(test7_indices);
R_evaluate_2(test7_indices) = R_rec_2(test7_indices);
R_origin (test7_indices) = R (test7_indices);

new_rating_0(row)=0;
new_rating(row)=0;
new_rating_2(row)=0;
ori_rating(row)=0;
precision_array_0(row)=0;
precision_array(row)=0;
precision_array_2(row)=0;
hit_array_0(row)=0;
hit_array(row)=0;
hit_array_2(row)=0;
alarm_array_0(row)=0;
alarm_array(row)=0;
alarm_array_2(row)=0;
j=1; 
Pre(1:20)=0;
Hit(1:20)=0;
Alarm(1:20)=0;
Pre_2(1:20)=0;
Hit_2(1:20)=0;
Alarm_2(1:20)=0;
Pre_0(1:20)=0;
Hit_0(1:20)=0;
Alarm_0(1:20)=0;


for L= 1:20
    for i= 1:row
    new_rating = R_evaluate(i,:);
    new_rating_2 = R_evaluate_2(i,:);
    new_rating_0 = R_evaluate_0(i,:);
    ori_rating = R_origin(i,:);
    [SortValue,SortIndex] = sort(new_rating(:),'descend');
    [SortValue_2,SortIndex_2] = sort(new_rating_2(:),'descend');
    [SortValue_0,SortIndex_0] = sort(new_rating_0(:),'descend');
    targetindex = SortIndex(1:L);
    targetindex_2 = SortIndex_2(1:L);
    targetindex_0 = SortIndex_0(1:L);
    precision_array(i)=length(find(ori_rating(targetindex)>3))/L;
    precision_array_2(i)=length(find(ori_rating(targetindex_2)>3))/L;
    precision_array_0(i)=length(find(ori_rating(targetindex_0)>3))/L;
    if length(find(ori_rating>3))>0
        hit_array(i)=length(find(ori_rating(targetindex)>3))/length(find(ori_rating>3));
        hit_array_2(i)=length(find(ori_rating(targetindex_2)>3))/length(find(ori_rating>3));
        hit_array_0(i)=length(find(ori_rating(targetindex_0)>3))/length(find(ori_rating>3));
    else
        hit_array(i)=0;
        hit_array_2(i)=0;
        hit_array_0(i)=0;
    end
    if length(find(ori_rating<3))>0
        alarm_array(i)=length(find(ori_rating(targetindex)<3))/length(find(ori_rating<3));
        alarm_array_2(i)=length(find(ori_rating(targetindex_2)<3))/length(find(ori_rating<3));
        alarm_array_0(i)=length(find(ori_rating(targetindex_0)<3))/length(find(ori_rating<3));
    else
        alarm_array(i)=0;
        alarm_array_2(i)=0;
        alarm_array_0(i)=0;
    end
    end

    Pre(j)=mean(precision_array);
    Hit(j)=sum(hit_array)/length(find(hit_array~=0));
    Alarm(j)=sum(alarm_array)/length(find(alarm_array~=0));
    Pre_2(j)=mean(precision_array_2);
    Hit_2(j)=sum(hit_array_2)/length(find(hit_array_2~=0));
    Alarm_2(j)=sum(alarm_array_2)/length(find(alarm_array_2~=0));
    Pre_0(j)=mean(precision_array_0);
    Hit_0(j)=sum(hit_array_0)/length(find(hit_array_0~=0));
    Alarm_0(j)=sum(alarm_array_0)/length(find(alarm_array_0~=0));
    j=j+1;
end
axisL=1:20;
subplot(3,1,1);
plot(axisL,Pre);
title('.*Rating version: Precision Rate');
subplot(3,1,2);
plot(axisL,Hit);
title('.*Rating version: Hit Rate');
subplot(3,1,3);
plot(axisL,Alarm);
title('.*Rating version: False Alarm Rate');

figure;
axisL=1:20;
subplot(3,1,1);
plot(axisL,Pre_2);
title('Regularized Precision Rate');
subplot(3,1,2);
plot(axisL,Hit_2);
title('Regularized Hit Rate');
subplot(3,1,3);
plot(axisL,Alarm_2);   
title('Regularized False Alarm Rate');
    
figure;
axisL=1:20;
subplot(3,1,1);
plot(axisL,Pre_0);
title('Basic version Precision Rate');
subplot(3,1,2);
plot(axisL,Hit_0);
title('Basic version Hit Rate');
subplot(3,1,3);
plot(axisL,Alarm_0);   
title('Basic version False Alarm Rate');

figure;
plot(Alarm,Hit,'r');hold on;
plot(Alarm_2,Hit_2,'b'); hold on;
plot(Alarm_0,Hit_0,'g');
title('ROC curve');
ylabel('hit rate');
xlabel('false alarm rate');
legend('Modified Method 2 with .*R','Regularized ALS','Method 2');
    
    