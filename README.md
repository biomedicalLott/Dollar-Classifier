# Dollar-Classifier
Designed by Robert Lott for Intelligent CAD - 
<p><img src = "https://i.imgur.com/LqCuRsD.png"/></p>
My implementation of a 1$ Recognizer for User Interface Prototypes

My dollar classifier, uses a simple template matching to distinguish between different shapes. I have implemented this within Matlab App which unfortunately at this time is quite slow. The processing is fast and I have written what I believe to be a fairly efficient implementation of the 1$ classifier. The problem is that matlab app doesn't take user input as quickly as I need it to for this classifier to work very well. At present this template matching's success will change pretty dramatically if you move from one computer to another with the app. 

With this app, I have built in the ability to create new templates from your drawings as well as to test the classifier. 

To create the template, the drawing you make is corrected for rotation and scaled to a very specific size. You can view the templates you create on the right hand side of the app. 
<p><img src = "https://i.imgur.com/Ake2hLp.png"/> </p>

Classification is based off of a score, the score is not normalized on output and as such will be some value positive or negative. 

Simply draw the shape in the space provided on the left 
<p><img src = "https://i.imgur.com/ICyFHnz.png"/> </p>

Then select "Classify"

Select clear between interactions so your gesture does not add to the data of the next interaction. 
