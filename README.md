# ML For Glucose Predicion

## Preparation

Please download the cleaned dataset using this [link](https://drive.google.com/drive/folders/1gLgYf21t9jJUanGhftExf9OHYGjVXy96?usp=sharing) and put it in the root directory. See below:

![Alt text](structure.png)
## Demo
one_module_demo.ipynb: an example of using KNN regression to predict the 10-week glucose level of mice, for a specific module in a particular tissue

demo.py: an example of using KNN regression to predict the 10-week glucose level of mice, FOR ALL MODULES in ALL TISSUES. A sample output file is presented in result/knn.csv

## Notes

- The number of mice across different tissues are not consistent

- We use 50 mice for testing, these 50 mice are the same across different tissues

- All the methods and notebooks about data cleaning are in notebooks/ ot utils.py.
