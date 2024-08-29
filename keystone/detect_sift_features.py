# Author: Shivesh Prakash
# This file contains code to detect SIFT feature points given two spectral images.

import numpy as np
import cv2


def detect_sift_features(image1, image2):
    """
    Detects SIFT (Scale-Invariant Feature Transform) features and matches them between two spectral images.

    Parameters:
    image1 (numpy.ndarray): The first spectral image, represented as a 2D numpy array.
    image2 (numpy.ndarray): The second spectral image, represented as a 2D numpy array.

    Returns:
    points1 (numpy.ndarray): Array of matched keypoint locations in the first image.
    points2 (numpy.ndarray): Array of matched keypoint locations in the second image.
    """

    # Normalize the spectral images to 8-bit grayscale.
    # SIFT expects input images in 8-bit format, so we scale the spectral values to fit this range.
    def normalize_to_uint8(image):
        normalized_img = cv2.normalize(image, None, 0, 255, cv2.NORM_MINMAX)
        return normalized_img.astype(np.uint8)

    # Apply normalization to both input images
    img1_normalized = normalize_to_uint8(image1)
    img2_normalized = normalize_to_uint8(image2)

    # Initialize SIFT detector using OpenCV's built-in function.
    sift = cv2.SIFT_create()

    # Detect keypoints and compute descriptors for both normalized images.
    keypoints1, descriptors1 = sift.detectAndCompute(img1_normalized, None)
    keypoints2, descriptors2 = sift.detectAndCompute(img2_normalized, None)

    # Match the descriptors between the two images using the FLANN (Fast Library for Approximate Nearest Neighbors) matcher.
    # FLANN is efficient for matching large descriptor sets.
    index_params = dict(
        algorithm=1, trees=5
    )  # Using KD-tree algorithm (algorithm=1) with 5 trees.
    search_params = dict(
        checks=50
    )  # Number of times the trees in the index should be recursively traversed.
    flann = cv2.FlannBasedMatcher(index_params, search_params)
    matches = flann.knnMatch(
        descriptors1, descriptors2, k=2
    )  # Find the 2 nearest neighbors for each descriptor.

    # Apply the ratio test to retain only good matches.
    # Lowe's ratio test is used to filter out false matches.
    good_matches = []
    for m, n in matches:
        if m.distance < 0.7 * n.distance:
            good_matches.append(m)

    # Extract the coordinates of the matched keypoints from both images.
    # These coordinates represent the locations of features that are common between the two images.
    points1 = np.float32([keypoints1[m.queryIdx].pt for m in good_matches])
    points2 = np.float32([keypoints2[m.trainIdx].pt for m in good_matches])

    return points1, points2
