// main.ts

import axios from "axios";

async function processFile() {
    const fileInput = document.getElementById('fileInput') as HTMLInputElement;
    const resultArea = document.getElementById('resultArea');

    const file = fileInput.files?.[0];

    if (file) {
        const formData = new FormData();
        formData.append('file', file);

        try {
            const response = await axios.post('/process', formData);
            const results = response.data;

            // Display results on the website
            resultArea.innerHTML = `<pre>${JSON.stringify(results, null, 2)}</pre>`;
        } catch (error) {
            console.error('Error processing file:', error);
        }
    } else {
        console.error('Please select a file.');
    }
}
