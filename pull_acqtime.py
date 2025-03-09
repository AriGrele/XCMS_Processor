import os
import xml.etree.ElementTree as ET
import csv

root_directory = 'F:/Dissertation chemistry/AG/Asclepias/Pollin/20240417'
output_csv = 'acqtime.csv'

with open(output_csv, 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(['Folder', 'AcqTime'])
    
    for i,folder_name in enumerate(os.listdir(root_directory)):
        print(i,folder_name)
        folder_path = os.path.join(root_directory, folder_name)
        
        if os.path.isdir(folder_path):
            xml_path = os.path.join(folder_path, 'AcqData', 'sample_info.xml')
            
            if os.path.exists(xml_path):
                try:
                    tree = ET.parse(xml_path)
                    root = tree.getroot()
                    
                    acq_time = None
                    for field in root.findall('.//Field'):
                        name = field.find('Name')
                        if name is not None and name.text == 'AcqTime':
                            value = field.find('Value')
                            if value is not None:
                                acq_time = value.text
                                break
                    
                    if acq_time:
                        csv_writer.writerow([folder_name, acq_time])
                except Exception as e:
                    print(f"Error processing {xml_path}: {e}")

print("done")
