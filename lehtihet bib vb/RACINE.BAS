Attribute VB_Name = "RACINE"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << RACINE >>
'  Requis   :
'  Externes : Fonc(nf,x)
'
'  Descrip. : Résolution d'une équation non linéaire
'
Option Explicit
Private Const mKMAX As Integer = 30
Function ESD(ByVal b As Double, _
             ByVal c As Double, _
             ByRef x1 As Double, _
             ByRef x2 As Double) As Integer
        
    'Résolution d'une équation du second degré x² + bx + c = 0
    ' La fonction détermine sigma le signe du discriminant
    ' Si sigma > 0 : x1 et x2 seront les deux solutions réelles
    '                et on aura toujours |x1| < |x2|
    ' Si Sigma = 0 : x1 sera la solution double
    '                x2 = 0 sera ignoré
    ' Si sigma < 0 : x1 et x2 seront respectivement
    '                les parties réelle et imaginaire
    '                des deux solutions complexes conjuguées
    
    Dim Sigma As Integer
    
    x1 = -0.5 * b
    x2 = x1 * x1 - c
    Sigma = Sgn(x2)
    x2 = Sqr(x2 * Sigma)
    If Sigma > 0 Then
        If x1 > 0 Then x2 = x1 + x2 Else x2 = x1 - x2
        x1 = c / x2
    End If
    ESD = Sigma
End Function
Function Dichotomie(ByVal nf As Integer, _
                    ByVal a As Double, _
                    ByVal b As Double, _
                    ByVal dx As Double, _
                    ByRef x As Double, _
           Optional ByRef y As Double) As Integer

    ' Placer ici la description + signification des paramètres
    
    Dim k As Integer, d As Double
    Dim xn As Double, xp As Double
    
    Call prv_cherche(nf, a, b, xn, xp)
    d = xp - xn ' taille de l'intervalle initial
    k = 0
    Do
        k = k + 1
        If k > mKMAX Then Dichotomie = -1: Exit Function ' non convergence
        d = 0.5 * d   ' Réduction de l'intervalle
        x = xp - d
        y = Fonc(nf, x)
        If y > 0 Then xp = x  ' Test d'encadrement
    Loop Until Abs(d) < dx Or y = 0   ' test d'arrêt
    Dichotomie = k + 2 ' nombre d'évaluation de la fonction
End Function
Private Sub prv_cherche(ByVal nf As Integer, _
                        ByVal a As Double, _
                        ByVal b As Double, _
                        ByRef xn As Double, _
                        ByRef xp As Double, _
               Optional ByRef yn As Double, _
               Optional ByRef yp As Double)
    
    ' Cette procédure cherche le côté positif
    ' et négatif d'une fonction externe Fonc
    ' définie par l'utilisateur.
    ' (xn,yn) et (xp,yp) sont les coordonnées
    ' des points extrêmes de l'intervalle [a,b]
    ' une erreur est signalée si Fonc(a) et Fonc(b)
    ' ne sont pas de signes opposés.
    
    yp = Fonc(nf, b)
    If yp > 0 Then
        xn = a
        xp = b
        yn = Fonc(nf, xn)
        If yn >= 0 Then Stop ' Données non valides
    ElseIf yp < 0 Then
        xn = b
        xp = a
        yn = yp
        yp = Fonc(nf, xp)
        If yp <= 0 Then Stop ' Données non valides
    Else
        Stop ' Données non valides
    End If
End Sub

 
