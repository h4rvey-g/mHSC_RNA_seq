#!/bin/bash

BAMLIST_DIR="data/102.splice_workflow/star_salmon/rmats/bamlist"
PREFIX="data/102.splice_workflow/star_salmon"

# 遍历所有bamlist文件
for bamlist in "$BAMLIST_DIR"/*_bamlist.txt; do
    # 获取输出文件名
    output_file="${bamlist%.*}_abs.txt"
    
    # 清空或创建输出文件
    > "$output_file"
    
    # 读取输入文件，添加前缀并写入新文件
    while IFS=, read -r line; do
        # 将行按逗号分割成数组
        IFS=',' read -ra files <<< "$line"
        
        # 初始化结果字符串
        result=""
        
        # 处理每个文件
        for ((i=0; i<${#files[@]}; i++)); do
            # 去除可能的空格
            file=${files[i]// /}
            
            # 添加前缀
            if [ $i -eq 0 ]; then
                result="$PREFIX/$file"
            else
                result="$result,$PREFIX/$file"
            fi
        done
        
        # 写入结果
        echo "$result" >> "$output_file"
    done < "$bamlist"
    
    echo "Created absolute path bamlist: $output_file"
done
