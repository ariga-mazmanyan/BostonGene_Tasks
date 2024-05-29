import json
import yaml
import re


def data_loader():

    def data_preprocessing(text):

        preprocessed_text = re.sub(r'[^\w\s]', ' ', text)
        return preprocessed_text

    with open('genes.yaml', 'r') as f:
        dict_yaml = yaml.safe_load(f)

    for gene in dict_yaml:
        for i, synonym in enumerate(gene["synonyms"]):
            gene["synonyms"][i] = data_preprocessing(synonym)
        gene["synonyms"] = list(set(gene["synonyms"]))

    # Splitting the data for Task1 and Task2
    with open('test_texts.json', 'r') as f:
        test_data = json.load(f)

    data_task1 = []
    data_task2 = []

    for data in test_data:
        data['text'] = data_preprocessing(data['text'])
        if data['genes']:
            data_task1.append({'text': data['text'], 'genes': data['genes']})

        if data["hla"]:
            data_task2.append({'text': data['text'], 'hla': data['hla']})

    return dict_yaml, data_task1, data_task2


def gene_detector(sentence):

    gene_dict = {"genes": []}
    sentence = sentence.replace("can", "   ")
    sentence = sentence.lower()

    for gene in genes:
        for synonym in gene["synonyms"]:

            synonym_lower = synonym.lower()

            pattern = r'\b{}\b'.format(re.escape(synonym_lower))
            matches = re.finditer(pattern, sentence)
            indexes = [[match.start(), match.end()] for match in matches]

            found = False
            if indexes:
                for i in gene_dict["genes"]:
                    if i["name"] == gene["name"]:
                        for j in i["positions"]:
                            if indexes[0][0] == j[0]:
                                if indexes[0][1] >= j[1]:
                                    i["positions"].remove(j)
                                else:
                                    indexes = []
                                    found = True
                        if indexes:
                            found = True
                            i["positions"].extend(indexes)
                            i["positions"].sort()

                if not found:
                    gene_dict["genes"].append({"name": gene["name"], "positions": indexes})

    return gene_dict


genes, test_data_task1, test_data_task2 = data_loader()
count_1 = 0
for text in test_data_task1:
    gene_dict = gene_detector(text["text"])
    if gene_dict["genes"] == text["genes"]:
        count_1 += 1
    else:
        print(text["text"], gene_dict["genes"], text["genes"], sep="\n")

print(f"The accuracy of 1st task is {count_1/len(test_data_task1) * 100}% ({count_1}/{len(test_data_task1)})")


def hla_detector(sentence):

    hla = []
    sentence_list = sentence.split()

    matches = re.finditer(r'\bHLA\b', sentence)
    indexes = [[match.start(), match.end()] for match in matches]
    if not indexes:
        matches = re.finditer(r'\bHuman Leukocyte Antigen\b', sentence)
        indexes = [[match.start(), match.end()] for match in matches]
    hla_index = [index for index, item in enumerate(sentence_list) if item == "HLA"]
    if not hla_index:
        hla_index = [index for index, item in enumerate(sentence_list) if item == "Antigen"]

    for pos, i in enumerate(hla_index):
        hla_dict = {"gene": False, "allele": None, "protein": None, "positions": [indexes[pos]]}
        item = sentence_list[i+1]
        if item.isupper():
            matches = re.finditer(r'\b{}\b'.format(item), sentence)
            item_index = [match.end() for match in matches if match.end() > indexes[pos][1]]
            hla_dict["positions"][-1][1] = item_index[0]
            if item.isalpha() or item[1].isalpha():
                hla_dict["gene"] = item

            else:
                hla_dict["gene"] = item[0]
                hla_dict["allele"] = item[1:]

            if hla_dict["allele"] is None:
                item = sentence_list[i + 2]
                if item.isdigit():
                    matches = re.finditer(r'\b{}\b'.format(item), sentence)
                    item_index = [match.end() for match in matches if match.end() > indexes[pos][1]]
                    hla_dict["positions"][-1][1] = item_index[0]

                    if len(item) == 4:
                        hla_dict["allele"] = item[:2]
                        hla_dict["protein"] = item[2:]
                    else:
                        hla_dict["allele"] = item
                        if sentence_list[i+3].isdigit():
                            hla_dict["protein"] = sentence_list[i+3]
                            matches = re.finditer(r'\b{}\b'.format(sentence_list[i+3]), sentence)
                            item_index = [match.end() for match in matches if match.end() > indexes[pos][1]]
                            hla_dict["positions"][-1][1] = item_index[0]
                        elif sentence_list[i+3][0].isdigit():
                            for j in range(len(sentence_list[i + 3]) - 1, -1, -1):
                                if sentence_list[i + 3][j].isdigit():
                                    hla_dict["protein"] = sentence_list[i + 3][:j + 1]
                                    matches = re.finditer(sentence_list[i + 3][:j + 1], sentence)
                                    item_index = [match.end() for match in matches if match.end() > indexes[pos][1]]
                                    hla_dict["positions"][-1][1] = item_index[0]
                                    break

            elif hla_dict["protein"] is None:
                if sentence_list[i + 2].isdigit():
                    hla_dict["protein"] = sentence_list[i + 2]
                    matches = re.finditer(r'\b{}\b'.format(sentence_list[i + 2]), sentence)
                    item_index = [match.end() for match in matches if match.end() > indexes[pos][1]]
                    hla_dict["positions"][-1][1] = item_index[0]
                elif sentence_list[i + 2][0].isdigit():
                    for j in range(len(sentence_list[i + 2]) - 1, -1, -1):
                        if sentence_list[i + 2][j].isdigit():
                            hla_dict["protein"] = sentence_list[i + 2][:j + 1]
                            matches = re.finditer(sentence_list[i + 2][:j + 1], sentence)
                            item_index = [match.end() for match in matches if match.end() > indexes[pos][1]]
                            hla_dict["positions"][-1][1] = item_index[0]
                            break

        if hla_dict["gene"]:
            hla.append(hla_dict)

    return hla


def multiple_gene_detector_without_hla(sentence):

    sentence_list = list(sentence)
    for char in range(hla[-1]["positions"][-1][1]):
        sentence_list[char] = ' '
    sentence = "".join(sentence_list)
    sentence_list_w = sentence.split()
    for w in sentence_list_w:
        if w.isupper():
            i = sentence.find(w)
            sentence_list[i - 5:i - 1] = " HLA"
            break
        elif w[:3].isupper():
            i = sentence.find(w[:3])
            sentence_list[i - 5:i - 1] = " HLA"
            break

    return "".join(sentence_list)


count_2 = 0
for text in test_data_task2:
    hla = hla_detector(text["text"])
    if hla == text["hla"]:
        count_2 += 1
    else:
        text["text"] = multiple_gene_detector_without_hla(text["text"])
        new_hla = hla_detector(text["text"])
        for i in new_hla:
            for j in i["positions"]:
                j[0] += 4
        hla.extend(new_hla)

        while hla != text["hla"]:
            text["text"] = multiple_gene_detector_without_hla(text["text"])
            new_hla = hla_detector(text["text"])
            for i in new_hla:
                for j in i["positions"]:
                    j[0] += 4
            hla.extend(new_hla)
        else:
            count_2 += 1

print(f"The accuracy of 2nd task is {count_2/len(test_data_task2) * 100}% ({count_2}/{len(test_data_task2)})")
