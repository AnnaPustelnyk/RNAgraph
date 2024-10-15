from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time

# Zainicjuj ChromeDriver bez konieczności podawania ścieżki
driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()))

try:
    # 1. Otwórz stronę
    driver.get("https://mcq.cs.put.poznan.pl/")

    # 2. Poczekaj na załadowanie strony
    time.sleep(3)

    # 3. Znajdź i kliknij przycisk wyboru "From local file"
    local_file_button = driver.find_element(By.XPATH, "//app-radio-button[@label='From local file']")
    local_file_button.click()

    # 4. Znajdź element do przesyłania pliku i załaduj plik
    upload_input = driver.find_element(By.XPATH, "//input[@type='file']")
    upload_input.send_keys("C:/Users/anna/Downloads/hello.pdb")  # Zmień na odpowiednią ścieżkę do twojego pliku

    # 5. Poczekaj na przetworzenie pliku
    time.sleep(5)

    # 6. Kliknij przycisk "submit", aby przetworzyć plik i sprawdzić interakcje
    submit_button = driver.find_element(By.XPATH, "//button[@color='primary']")
    submit_button.click()

    # 7. Poczekaj na wyniki (aż tabela z interakcjami się pojawi)
    WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.CLASS_NAME, 'data-panel'))
    )

    # 8. Wyciągnij informacje o interakcjach
    rows = driver.find_elements(By.XPATH, "//mat-row[contains(@class, 'cdk-row')]")  # Wszystkie wiersze tabeli

    for row in rows:
        # Pobierz informacje o nukleotydach i typie interakcji
        nucleotides = row.find_element(By.XPATH, ".//mat-cell[contains(@class, 'cdk-column-basePair')]").text
        interaction_type = row.find_element(By.XPATH, ".//mat-cell[contains(@class, 'cdk-column-interactionType')]").text
        saenger_code = row.find_element(By.XPATH, ".//mat-cell[contains(@class, 'cdk-column-saenger')]").text
        
        # Wyświetl uzyskane informacje
        print(f"Nukleotydy: {nucleotides}, Typ interakcji: {interaction_type}, Kod Saengera: {saenger_code}")

finally:
    # 9. Zamknij przeglądarkę
    driver.quit()
